#!/usr/bin/env python3

"""\
Plot which residues were chosen to be deleted.

Usage:
    minp plot_deletions <dels_workspace>... [-n]

Options:
    -n --normalize
        Scale each plot by the total number of deletions, to make it easier to 
        compare strategies that produce very different numbers of deletions.  
        Note that strategies that produce longer deletions will still have more 
        area than those that produce shorter deletions.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from minp import DeletionsWorkspace, load_weighted_msa

def main():
    import docopt
    args = docopt.docopt(__doc__)

    for path in args['<dels_workspace>']:
        work_dels = DeletionsWorkspace.from_path(path)
        msa = load_weighted_msa(work_dels.msa)
        dels = pd.read_hdf(work_dels.deletions_hdf5)

        y = count_deletions(msa, dels)
        x = np.arange(len(y))
        label = work_dels.root.relative_to(work_dels.msa.blast.root.parent)

        plt.plot(x, y, label=f'{label} (N={len(dels)})')

    plt.xlabel("residue index")
    plt.ylabel("relative deletions" if args['--normalize'] else "deletions")
    plt.legend(loc='best')
    plt.show()

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


