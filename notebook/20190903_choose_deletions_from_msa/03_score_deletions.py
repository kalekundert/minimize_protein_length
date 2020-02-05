#!/usr/bin/env python3

"""\
Usage:
    03_calc_deletion_scores <msa_workspace>
"""


import docopt
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from Bio.SubsMat import MatrixInfo
from workspace import MsaWorkspace, load_msa, weight_alignments

def calc_deletion_scores(ref, msa):
    """
    Sum the weight of each deletion at each position in the reference 
    sequence.

    My goal is for this score to be independent of the number of sequences in 
    the alignment.  Adding good alignments shouldn't change things down-weight 
    less common deletions too much, and adding poor alignments should have a 
    negligible effect on the scores.
    """
    i = map_msa_indices_to_ref(ref, msa)
    deletion_scores = np.zeros(len(i))
    all_weights = 0

    for record in msa:
        for msa_i in i:
            all_weights += record.weight
            deletion_scores[i[msa_i]] += \
                    record.weight * (record.seq[msa_i] == '-')

    return deletion_scores / all_weights

def calc_conservation_scores(ref, msa):
    i = map_msa_indices_to_ref(ref, msa)
    conservation_scores = np.zeros(len(i))
    subs_mat = MatrixInfo.blosum80
    all_weights = 0

    def apply_subs_mat(a, b):
        if (a,b) in subs_mat:
            return subs_mat[a,b]
        if (b,a) in subs_mat:
            return subs_mat[b,a]
        return 0

    for msa_i in i:
        for record in msa:
            subs_score = apply_subs_mat(
                    ref.seq[msa_i],
                    record.seq[msa_i],
            )

            all_weights += record.weight
            conservation_scores[i[msa_i]] += record.weight * subs_score

    return conservation_scores / all_weights

def map_msa_indices_to_ref(ref, msa):
    ref_i_from_msa_j = {}

    ref_i = 0
    for msa_j, seq in enumerate(ref.seq):
        if seq == '-': continue
        ref_i_from_msa_j[msa_j] = ref_i
        ref_i += 1

    return ref_i_from_msa_j


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    work_blast, work_msa = MsaWorkspace.from_path(args['<msa_workspace>'])

    # Weight the alignment and report statistics.
    ref, msa = load_msa(work_msa.output_aln, work_blast.query)
    weight_alignments(ref, msa)

    print(f"{len(msa)} sequences aligned")
    print(f"{sum(x.percent_id > 0.3 for x in msa)} exceed 30% identity")

    k = gaussian_kde([x.percent_id for x in msa])
    x = np.linspace(0, 1, num=200)

    plt.subplot(2, 1, 1)
    plt.plot(x, k(x))
    plt.axvline(0.3, color='red')
    plt.ylabel('Count')
    plt.xlabel('% Identity')
    
    # Calculate and record the deletion scores.
    deletion_scores = calc_deletion_scores(ref, msa)
    deletion_scores.tofile(work_msa.deletion_scores_cache)

    threshold = np.mean(deletion_scores)
    n = len(deletion_scores)
    n_del = sum(deletion_scores > threshold)
    print(f"{n_del}/{n} residues ({100*n_del/n:.1f}%) have above-average deletion scores.")

    i = np.arange(n)
    plt.subplot(2, 1, 2)
    plt.plot(i, deletion_scores)
    plt.tight_layout()
    plt.ylabel('Score')
    plt.xlabel('Residue index')
    plt.xlim((0, len(i)))
    plt.axhline(threshold, color='red')

    plt.title(work_msa.root.name)
    plt.savefig(work_msa.deletion_scores_plot)
    plt.show()
