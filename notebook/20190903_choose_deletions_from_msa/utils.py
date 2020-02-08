#!/usr/bin/env python3

import shutil, time, re
import autoprop, toml
import numpy as np
from pathlib import Path
from datetime import date, timedelta
from functools import lru_cache
from inform import fatal
from copy import deepcopy
from contextlib import contextmanager

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.SubsMat import MatrixInfo

@autoprop
class Workspace:

    def __init__(self, root):
        self.root = Path(root).resolve()
        self.metadata_path = self.root / 'metadata.toml'
        self.mkdir()

    def mkdir(self):
        self.root.mkdir(parents=True, exist_ok=True)

    def rmdir(self):
        if self.root.exists():
            shutil.rmtree(self.root)

    @contextmanager
    def touch(self, path):
        yield path
        self.metadata['modification_dates'][path.name] = date.today()

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
            fatal("{self.repr} parameters differ from those used previously!")
            fatal()
            for key in params:
                prev = prev_params[key]
                curr = curr_params[key]
                if prev != curr:
                    codicil(f"    {key} was {prev!r}, now {curr!r}")
            fatal()
            fatal("Use the -f flag to overwrite.  Aborting.")
            sys.exit(1)

        self.metadata['parameters'] = curr_params

@autoprop
class BlastWorkspace(Workspace):
    repr = 'BLAST'

    def __init__(self, root):
        super().__init__(root)
        self.output_xml = self.root / 'blast_hits.xml'
        self.output_fasta = self.root / 'blast_hits.fasta'

    def get_query(self):
        return self.metadata['parameters']['query']

    def get_database(self):
        return self.metadata['parameters']['database']

    def get_count(self):
        return self.metadata['results']['num_hits']

class MsaWorkspace(Workspace):
    repr = 'MSA'

    @classmethod
    def from_params(cls, work_blast, algorithm, limit):
        name = f'{algorithm}_{limit}' if limit else algorithm
        return cls(work_blast, name)

    @classmethod
    def from_path(cls, root):
        root = Path(root)
        work_blast = BlastWorkspace(root.parent)
        return work_blast, cls(work_blast, root.name)


    def __init__(self, work_blast, name):
        super().__init__(work_blast.root / name)
        self.blast = work_blast

        self.input_fasta = self.root / 'input.fasta'
        self.output_aln = self.root / 'output.aln'
        self.pruned_aln = self.root / 'output_pruned.aln'
        self.deletion_scores_cache = self.root / 'deletion_scores.npy'
        self.deletion_scores_plot = self.root / 'deletion_scores.svg'

class DeletionsWorkspace(Workspace):

    @classmethod
    def from_path(cls, root):
        root = Path(root)
        work_blast, work_msa = MsaWorkspace.from_path(root.parent)
        return cls(work_msa, root.name)


    def __init__(self, work_msa, name):
        super().__init__(work_msa.root / name)
        self.msa = work_msa

        self.deletions_fasta = self.root / 'deletions.fasta'
        self.deletions_aln = self.root / 'deletions.aln'
        self.deletions_hdf5 = self.root / 'deletions.hdf5'

    def write_deletions(self, dels, msa):
        """
        dels is a data frame with 'del_start' and 'del_end' columns.
        """
        records = []
        aligned_records = []

        for i, row in dels.iterrows():
            # Because all the columns in the dataframe are numeric, pandas 
            # automatically coerces the ints to floats.  I think this is a bug (see 
            # pandas #28308), but it's easy to work around.
            start, end = int(row['del_start']), int(row['del_end'])
            start1 = start + 1  # 1-indexed for human-readable output
            gap = end - start
            assert end <= len(msa.ref_ungapped), row

            id = f'{msa.ref.id}Δ{start1}' + (f'-{end}' if start1 != end else '')
            desc = ' '.join(
                    f'{k}={repr(v)}'
                    for k, v in row.items()
                    if k not in {'del_start', 'del_end'}
            )
            seq = Seq(
                msa.ref_ungapped[:start] + msa.ref_ungapped[end:],
                IUPAC.protein,
            )
            aligned_seq = Seq(
                msa.ref_ungapped[:start] + ('-'*gap) + msa.ref_ungapped[end:],
                Gapped(IUPAC.protein, '-'),
            )

            record = SeqRecord(seq, id=id, description=desc)
            records.append(record)

            aligned_record = SeqRecord(aligned_seq, id=id, description=desc)
            aligned_records.append(aligned_record)

        # FASTA file to design oligos.
        with self.touch(self.deletions_fasta) as p:
            SeqIO.write(records, p, "fasta")

        # MSA to view in snapgene.
        msa = MultipleSeqAlignment(aligned_records)
        with self.touch(self.deletions_aln) as p:
            AlignIO.write(msa, p, 'clustal')

        # HDF5 to load in subsequent analysis scripts.
        with self.touch(self.deletions_hdf5) as p:
            dels.to_hdf(p, 'dels')

class LoophashWorkspace(DeletionsWorkspace):

    def __init__(self, work_msa):
        super().__init__(work_msa, 'del_loophash')

        self.input_pdb = self.msa.blast.root / '4un3.cif'
        self.input_pdb_chain = 'B'
        self.input_loophash_db = self.msa.blast.root / 'loophash_db'

        self.recovery_hdf5 = self.root / 'gap_recovery.hdf5'
        self.spannable_hdf5 = self.root / 'spannable_gaps.hdf5'
        self.filters_toml = self.root / 'spannable_gap_filters.toml'
        self.rosetta_log = self.root / 'rosetta_log.txt'

def load_weighted_msa(work_msa):
    """
    The given multiple sequence alignment (MSA) should contain the reference 
    sequence.  The reference will be removed from the alignment and returned 
    separately.  The alignment will also be changed such that "." is used to 
    indicate terminal deletions while "-" is used to indicate internal 
    deletions.
    """
    msa_with_ref = AlignIO.read(work_msa.output_aln, 'clustal')

    ref = None
    msa = MultipleSeqAlignment([])

    for record in msa_with_ref:
        # Use "." to indicate terminal mismatches, and "-" to indicate internal 
        # mismatches.

        to_dots = lambda m: '.' * (m.end() - m.start())
        record.seq = Seq(
                re.sub('^-+|-+$', to_dots, str(record.seq)),
                record.seq.alphabet,
        )

        if record.id == work_msa.blast.query:
            ref = record
        else:
            msa.append(record)

    msa.ref = ref
    msa.ref_ungapped = remove_gaps(ref.seq)

    weight_alignments(msa)

    return msa

def weight_alignments(msa):
    """
    Weight each alignment by percent identity to a reference sequence.

    The reference sequence is expected to be found in the 'ref' attribute of 
    the given MSA.  This sequence should be aligned to the MSA (i.e. contain 
    gaps), such that each index in the sequence is comparable to that same 
    index in the MSA.
    """

    def percent_identity(a, b):
        assert len(a) == len(b)
        N = len(b.strip('.'))
        return sum(aa == bb and aa not in '-.' for aa, bb in zip(a, b)) / N

    def scale(x, cutoff=0.3):
        m = 1 / (1 - cutoff)
        b = - m * cutoff
        return max(m * x + b, 0)

    for record in msa:
        record.percent_id = percent_identity(msa.ref.seq, record.seq)

        # I decided not to scale the weights, because it reduces the diversity 
        # of the results.
        record.weight = record.percent_id

    msa.sort(key=lambda x: x.weight, reverse=True)

def calc_deletion_scores(msa):
    """
    Sum the weight of each deletion at each position in the reference 
    sequence.

    My goal is for this score to be independent of the number of sequences in 
    the alignment.  Adding good alignments shouldn't change things down-weight 
    less common deletions too much, and adding poor alignments should have a 
    negligible effect on the scores.
    """

    i = map_ungapped_indices(msa.ref.seq)
    deletion_scores = np.zeros(len(i))

    for record in msa:
        for msa_i in i:
            deletion_scores[i[msa_i]] += \
                    record.weight * (record.seq[msa_i] == '-')

    return deletion_scores

def calc_conservation_scores(msa):
    i = map_ungapped_indices(msa.ref.seq)
    conservation_scores = np.zeros(len(i))
    subs_mat = MatrixInfo.blosum80

    def apply_subs_mat(a, b):
        if (a,b) in subs_mat:
            return subs_mat[a,b]
        if (b,a) in subs_mat:
            return subs_mat[b,a]
        return 0

    for msa_i in i:
        for record in msa:
            subs_score = apply_subs_mat(
                    msa.ref.seq[msa_i],
                    record.seq[msa_i],
            )

            conservation_scores[i[msa_i]] += record.weight * subs_score

    return conservation_scores

def map_ungapped_indices(seq):
    ungapped_indices = {}

    i = 0
    for j, aa in enumerate(seq):
        if aa in '-.': continue
        ungapped_indices[j] = i
        i += 1

    return ungapped_indices

def remove_gaps(seq):
    unalign = str.maketrans('', '', '-.')
    return str(seq).translate(unalign)

@contextmanager
def report_elapsed_time():
    start = time.time()
    yield
    end = time.time()
    print(f"Finished in {timedelta(seconds=end - start)}\n")

import pytest

@pytest.mark.parametrize(
        'seq, expected', [(
            'ABC',
            'ABC',
        ), (
            '.BC',
            'BC',
        ), (
            'A-C',
            'AC',
        ), (
            'AB.',
            'AB',
        )]
)
def test_remove_gaps(seq, expected):
    assert remove_gaps(seq) == expected

@pytest.mark.parametrize(
        'seq, expected', [(
            'AAA',
            {0:0, 1:1, 2:2},
        ), (
            '.AA',
            {1:0, 2:1},
        ), (
            'A-A',
            {0:0, 2:1},
        ), (
            'AA.',
            {0:0, 1:1},
        )]
)
def test_ungapped_indices(seq, expected):
    assert map_ungapped_indices(seq) == expected

