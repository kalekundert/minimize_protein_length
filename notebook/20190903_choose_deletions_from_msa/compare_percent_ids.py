#!/usr/bin/env python3

"""\
Compare distributions of percent identity to the reference sequence for 
multiple MSAs.

Usage:
    03_compare_percent_ids.py <msa_workspaces>...
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from utils import MsaWorkspace, load_msa, weight_alignments

if __name__ == '__main__':

    # Parse the command-line arguments.

    import docopt
    args = docopt.docopt(__doc__)

    # Make a histogram for each given MSA.

    for path in args['<msa_workspaces>']:
        work_blast, work_msa = MsaWorkspace.from_path(path)
        ref, msa = load_msa(work_msa.output_aln, work_blast.query)
        weight_alignments(ref, msa)


        k = gaussian_kde([x.percent_id for x in msa])
        x = np.linspace(0, 1, num=200)

        plt.plot(x, k(x), label=work_msa.root.name)

    # Show the histograms.

    plt.ylabel('Count')
    plt.xlabel('% Identity')
    plt.axvline(0.3, color='grey', linestyle='--', zorder=-1)
    plt.legend(loc='best')
    plt.show()
