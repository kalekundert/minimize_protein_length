#!/usr/bin/env python3

import shutil, time, re
import autoprop, toml
from pathlib import Path
from datetime import date, timedelta
from functools import lru_cache
from inform import fatal
from copy import deepcopy
from contextlib import contextmanager

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
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

def load_msa(msa_path, ref_id):
    """
    The given multiple sequence alignment (MSA) should contain the reference 
    sequence.  The reference will be removed from the alignment and returned 
    separately.  The alignment will also be changed such that "." is used to 
    indicate terminal deletions while "-" is used to indicate internal 
    deletions.
    """
    msa_with_ref = AlignIO.read(msa_path, 'clustal')

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

        if record.id == ref_id:
            ref = record
        else:
            msa.append(record)

    return ref, msa

def weight_alignments(ref, msa):
    """
    Weight each alignment by percent identity to the reference.

    The weight is scaled such that anything below 30% identity has zero weight, 
    because 30% is the level at which automatic MSA algorithms are considered 
    to perform poorly.
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
        record.percent_id = percent_identity(ref.seq, record.seq)

        # I decided not to scale the weights, because it reduces the diversity 
        # of the results.
        record.weight = record.percent_id

    msa.sort(key=lambda x: x.weight, reverse=True)

@contextmanager
def report_elapsed_time():
    start = time.time()
    yield
    end = time.time()
    print(f"Finished in {timedelta(seconds=end - start)}\n")
