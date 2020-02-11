#!/usr/bin/env python3

import shutil, time
import autoprop, toml
from pathlib import Path
from datetime import date, timedelta
from inform import display, fatal
from contextlib import contextmanager
from more_itertools import one

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped, generic_dna

class Workspace:

    def __init__(self, root):
        self.root = Path(root).resolve()

    def mkdir(self):
        self.root.mkdir(parents=True, exist_ok=True)

    def rmdir(self):
        if self.root.exists():
            shutil.rmtree(self.root)

    def exists(self):
        return self.root.exists()

@autoprop
class SharedWorkspace(Workspace):

    def __init__(self, root):
        super().__init__(root)
        self.params_path = self.root / 'params.toml'
        self.target_fasta = self.root / 'target.fasta'
        self.loophash_db = self.root / 'loophash_db'

    @autoprop.cache
    def get_params(self):
        if self.params_path.exists():
            return toml.load(self.params_path)
        else:
            return dict(
                    target={},
                    user={},
            )

    def get_target_id(self):
        return self.target_dna_record.id

    @autoprop.cache
    def get_target_dna_record(self):
        return one(SeqIO.parse(self.target_fasta, 'fasta'))

    def get_target_dna_seq(self):
        return self.target_dna_record.seq

    @autoprop.cache
    def get_target_protein_record(self):
        dna = self.target_dna_record
        return SeqRecord(
                dna.seq.translate(),
                id=dna.id,
                name=dna.name,
                description=dna.description,
        )

    def get_target_protein_seq(self):
        return self.target_protein_record.seq

    def get_target_pdb(self, src=None):
        if src:
            suffix = Path(src).suffix
            return self.root / f'target{suffix}'
        else:
            for suffix in ['.pdb', '.cif']:
                p = self.root / f'target{suffix}'
                if p.exists():
                    return p
            raise ValueError

    def get_target_pdb_chain(self):
        return self.params['target']['pdb_chain']

    def get_user_email(self):
        return self.params['user']['email']

    def write_params(self):
        with self.params_path.open('w') as f:
            toml.dump(self.params, f)

@autoprop
class JobWorkspace(Workspace):

    def __init__(self, work_shared, root):
        super().__init__(root)
        self.shared = work_shared
        self.metadata_path = self.root / 'metadata.toml'
        self.mkdir()

    @contextmanager
    def touch(self, path):
        yield path
        self.metadata['modification_dates'][path.name] = date.today()

    def get_relpath(self):
        return self.root.relative_to(self.shared.root.parent)

    @autoprop.cache
    def get_metadata(self):
        if self.metadata_path.exists():
            return toml.load(self.metadata_path)
        else:
            return dict(
                    parameters={},
                    modification_dates={},
            )

    def write_metadata(self):
        with self.metadata_path.open('w') as f:
            toml.dump(self.metadata, f)

    def update_params(self, **params):
        prev_params = self.metadata['parameters']
        curr_params = params

        if prev_params and prev_params != curr_params:
            error("{self.repr} parameters differ from those used previously!")
            error()
            for key in params:
                prev = prev_params[key]
                curr = curr_params[key]
                if prev != curr:
                    codicil(f"    {key} was {prev!r}, now {curr!r}")
            error()
            error("Use the -f flag to overwrite.  Aborting.")
            sys.exit(1)

        self.metadata['parameters'] = curr_params

@autoprop
class HomologsWorkspace(JobWorkspace):
    repr = 'BLAST'

    @classmethod
    def from_params(cls, shared_path, *, algorithm, database, num_hits):
        work_in = SharedWorkspace(shared_path)
        database_slug = database.split('_')[0]
        name = f'{algorithm}_{database_slug}_{num_hits}'
        return cls(work_in, name)

    @classmethod
    def from_path(cls, path):
        root = Path(path).resolve()
        work_in = SharedWorkspace(root.parent)
        return cls(work_in, root.name)


    def __init__(self, work_in, name):
        super().__init__(work_in, work_in.root / name)
        self.input_fasta = self.root / 'query.fasta'
        self.output_xml = self.root / 'homologs.xml'
        self.output_fasta = self.root / 'homologs.fasta'

    def get_algorithm(self):
        return self.metadata['parameters']['algorithm']

    def get_query(self):
        return str(self.shared.target_protein_seq)

    def get_database(self):
        return self.metadata['parameters']['database']

    def get_num_hits(self):
        return self.metadata['parameters']['num_hits']

    @autoprop.cache
    def get_output_records(self):
        return list(SeqIO.parse(self.output_fasta, 'fasta'))

@autoprop
class MsaWorkspace(JobWorkspace):
    repr = 'MSA'

    @classmethod
    def from_params(cls, homologs_path, *, algorithm, limit):
        work_homs = HomologsWorkspace.from_path(homologs_path)
        name = f'{algorithm}_{limit}' if limit else algorithm
        return cls(work_homs, name)

    @classmethod
    def from_path(cls, root):
        root = Path(root).resolve()
        work_homs = HomologsWorkspace.from_path(root.parent)
        return work_homs, cls(work_homs, root.name)


    def __init__(self, work_homs, name):
        super().__init__(work_homs.shared, work_homs.root / name)
        self.homologs = work_homs

        self.input_fasta = self.root / 'input.fasta'
        self.output_aln = self.root / 'output.aln'
        self.pruned_aln = self.root / 'output_pruned.aln'

    def get_algorithm(self):
        return self.metadata['parameters']['algorithm']

    def get_limit(self):
        return self.metadata['parameters']['limit']

class DeletionsWorkspace(JobWorkspace):

    @classmethod
    def from_path(cls, root):
        root = Path(root).resolve()
        work_homs, work_msa = MsaWorkspace.from_path(root.parent)
        return cls(work_msa, root.name)


    def __init__(self, work_msa, name):
        super().__init__(work_msa.shared, work_msa.root / name)
        self.msa = work_msa

        self.deletions_fasta = self.root / 'deletions.fasta'
        self.deletions_aln = self.root / 'deletions.aln'
        self.deletions_hdf5 = self.root / 'deletions.hdf5'

    def write_deletions(self, dels):
        """
        dels is a data frame with 'del_start' and 'del_end' columns.
        """
        deletion_records = []
        aligned_records = []

        ref_dna = self.shared.target_dna_seq
        ref_protein = self.shared.target_protein_seq

        for i, row in dels.iterrows():
            # Because all the columns in the dataframe are numeric, pandas 
            # automatically coerces the ints to floats.  I think this is a bug 
            # (see pandas #28308), but it's easy to work around.
            start, end = int(row['del_start']), int(row['del_end'])
            start1 = start + 1  # 1-indexed for human-readable output
            gap = end - start
            assert end <= len(ref_protein), row

            id = f'{self.shared.target_id}Δ{start1}' + (f'-{end}' if start1 != end else '')
            desc = f"score={row['del_score']:.1f}"
            deletion_seq = ref_dna[3*start:3*end]
            aligned_seq = Seq(
                str(ref_protein[:start] + ('-' * gap) + ref_protein[end:]),
                Gapped(IUPAC.protein, '-'),
            )

            deletion_record = SeqRecord(deletion_seq, id=id, description=desc)
            deletion_records.append(deletion_record)

            aligned_record = SeqRecord(aligned_seq, id=id, description=desc)
            aligned_records.append(aligned_record)

        # HDF5 to load in subsequent analysis scripts.
        with self.touch(self.deletions_hdf5) as p:
            dels.to_hdf(p, 'dels')

        # FASTA file to design oligos.
        with self.touch(self.deletions_fasta) as p:
            SeqIO.write(deletion_records, p, "fasta")

        # MSA to view in snapgene.
        msa = MultipleSeqAlignment(aligned_records)
        with self.touch(self.deletions_aln) as p:
            AlignIO.write(msa, p, 'clustal')

class LoophashWorkspace(DeletionsWorkspace):

    def __init__(self, work_msa):
        super().__init__(work_msa, 'loophash')

        self.recovery_hdf5 = self.root / 'gap_recovery.hdf5'
        self.spannable_hdf5 = self.root / 'spannable_gaps.hdf5'
        self.filters_toml = self.root / 'spannable_gap_filters.toml'
        self.rosetta_log = self.root / 'rosetta_log.txt'

@contextmanager
def report_elapsed_time():
    start = time.time()
    yield
    end = time.time()
    delta = timedelta(seconds=end - start)
    display(f"Finished in {delta}\n")

