#!/usr/bin/env python3

"""\
Compare distributions of percent identity to the reference sequence for 
multiple MSAs.

Usage:
    minp compare_percent_ids <msa_workspaces>...
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from minp import MsaWorkspace, load_weighted_msa, weight_alignments

def main():

    # Parse the command-line arguments.

    import docopt
    args = docopt.docopt(__doc__)

    # Make a histogram for each given MSA.

    work_msas = [
            MsaWorkspace.from_path(p)[1]
            for p in args['<msa_workspaces>']
    ]
    by_limit_then_algorithm = lambda x: (x.limit, x.algorithm)

    for work_msa in sorted(work_msas, key=by_limit_then_algorithm):
        msa = load_weighted_msa(work_msa)

        y, xx = np.histogram(
                [x.percent_id for x in msa],
                bins=100,
                range=(0,1),
        )
        x = (xx[:-1] + xx[1:]) / 2
        plt.plot(x, y, label=work_msa.relpath)


    # Show the histograms.

    plt.ylabel('Count')
    plt.xlabel('% Identity')
    plt.axvline(0.3, color='grey', linestyle='--', zorder=-1)
    plt.legend(loc='best')
    plt.title(work_msa.shared.root.name)
    plt.gcf().canvas.set_window_title(work_msa.shared.root.name)
    plt.show()
