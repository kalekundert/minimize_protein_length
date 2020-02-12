#!/usr/bin/env python3

import re
import numpy as np

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SubsMat import MatrixInfo

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

        if record.id == work_msa.shared.target_id:
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

    N = len(remove_gaps(msa.ref.seq))

    def percent_identity(a, b):
        assert len(a) == len(b)
        return sum(ai == bi and ai not in '-.' for ai, bi in zip(a, b)) / N

    for record in msa:
        assert len(msa.ref.seq) == len(record.seq)
        record.percent_id = percent_identity(msa.ref.seq, record.seq)

        # The weight is the same as the percent identity.  I made two variables 
        # because I was playing with some different weighting schemes in the 
        # past, but settled on simply having the two be the same.
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

def count_deletions(msa, dels, normalize=False):
    n_res = len(msa.ref_ungapped)
    n_dels = np.zeros(n_res, dtype=int)

    for _, row in dels.iterrows():
        i = int(row['del_start'])
        j = int(row['del_end'])
        n_dels[i:j] += 1

    if normalize:
        n_dels = n_dels / sum(n_dels)

    return n_dels

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

