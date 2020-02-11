#!/usr/bin/env python3

"""\
Plot the deletions scores for each residue in the target protein.

Usage:
    minp plot_deletion_scores <msa_workspace>
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from Bio.SubsMat import MatrixInfo
from minp import MsaWorkspace, load_weighted_msa, calc_deletion_scores

def main():

    # Parse the command-line arguments.

    import docopt
    args = docopt.docopt(__doc__)
    work_homs, work_msa = MsaWorkspace.from_path(args['<msa_workspace>'])

    msa = load_weighted_msa(work_msa)
    scores = calc_deletion_scores(msa)

    if os.fork():
        sys.exit()

    # Report statistics about the MSA.

    print(f"{len(msa)} sequences aligned")
    print(f"{sum(x.percent_id > 0.3 for x in msa)} exceed 30% identity")

    y, xx = np.histogram(
            [x.percent_id for x in msa],
            bins=100,
            range=(0,1),
    )
    x = (xx[:-1] + xx[1:]) / 2

    plt.subplot(2, 1, 1)
    plt.title(work_msa.relpath)
    plt.plot(x, y)
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

    plt.gcf().canvas.set_window_title(str(work_msa.relpath))
    plt.show()
