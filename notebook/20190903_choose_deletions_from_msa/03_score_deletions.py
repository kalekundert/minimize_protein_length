#!/usr/bin/env python3

"""\
Usage:
    03_calc_deletion_scores <msa_workspace>
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from Bio.SubsMat import MatrixInfo
from utils import *

if __name__ == '__main__':

    # Parse the command-line arguments.

    import docopt
    args = docopt.docopt(__doc__)
    work_blast, work_msa = MsaWorkspace.from_path(args['<msa_workspace>'])

    # Weight the alignment and report statistics.

    ref, msa = load_weighted_msa(work_msa)
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
