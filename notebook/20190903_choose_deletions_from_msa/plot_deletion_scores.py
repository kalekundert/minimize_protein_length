#!/usr/bin/env python3

"""\
Usage:
    plot_deletion_scores.py <msa_workspace>
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

    msa = load_weighted_msa(work_msa)
    scores = calc_deletion_scores(msa)

    # Report statistics about the MSA.

    print(f"{len(msa)} sequences aligned")
    print(f"{sum(x.percent_id > 0.3 for x in msa)} exceed 30% identity")

    k = gaussian_kde([x.percent_id for x in msa])
    x = np.linspace(0, 1, num=200)

    plt.subplot(2, 1, 1)
    plt.title(work_msa.root.relative_to(work_msa.blast.root.parent))
    plt.plot(x, k(x))
    plt.axvline(0.3, color='red')
    plt.ylabel('Count')
    plt.xlabel('% Identity')
    
    # Report statistics about the deletion scores.

    threshold = np.mean(scores)
    n = len(scores)
    n_del = sum(scores > threshold)
    print(f"{n_del}/{n} residues ({100*n_del/n:.1f}%) have above-average deletion scores.")

    i = np.arange(n)
    plt.subplot(2, 1, 2)
    plt.plot(i, scores)
    plt.tight_layout()
    plt.ylabel('Score')
    plt.xlabel('Residue index')
    plt.xlim((0, len(i)))
    plt.axhline(threshold, color='red')

    plt.show()
