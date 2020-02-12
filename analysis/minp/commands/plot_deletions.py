#!/usr/bin/env python3

"""\
Plot which residues were chosen to be deleted.

Usage:
    minp plot_deletions <dels_workspace>... [-n] [-a]

Options:
    -n --normalize
        Scale each plot by the total number of deletions, to make it easier to 
        compare strategies that produce very different numbers of deletions.  
        Note that strategies that produce longer deletions will still have more 
        area than those that produce shorter deletions.

    -a --align
        When showing scores for different sequences, align the sequences so 
        that aligned residues are plotted at the same position on the x-axis.  
        This requires that exactly two workspaces are provided.
"""

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from minp import DeletionsWorkspace, load_weighted_msa, count_deletions
from dataclasses import dataclass
from inform import fatal

from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat.MatrixInfo import blosum62

@dataclass
class Plot:
    x: np.array = None
    y: np.array = None
    label: str = ''
    seq: str = ''

def main():
    import docopt
    args = docopt.docopt(__doc__)

    plots = []

    for path in args['<dels_workspace>']:
        work_dels = DeletionsWorkspace.from_path(path)
        msa = load_weighted_msa(work_dels.msa)
        dels = pd.read_hdf(work_dels.deletions_hdf5)

        plot = Plot()
        plot.y = count_deletions(msa, dels)
        plot.x = np.arange(len(plot.y))
        plot.label = f'{work_dels.relpath} (N={len(dels)})'
        plot.seq = msa.ref_ungapped
        plots.append(plot)

    if args['--align']:
        if len(plots) != 2:
            fatal("Must specify 2 worksapces to use the --align option.")

        # I decided to use BLOSUM62 because the two sequences in this case may 
        # not be particularly similar.  I used the corresponding gap penalties 
        # from BLAST 2.2.27, which I found in the reference below:
        #
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/

        alignments = align.globalds(
                plots[0].seq,
                plots[1].seq,
                blosum62,
                -11, -1,
        )
        aligned_seq1, aligned_seq2, score, start, end = alignments[0]
        aligned_x1 = []
        aligned_x2 = []

        for i, (aa1, aa2) in enumerate(zip(aligned_seq1, aligned_seq2)):
            if aa1 not in '-':
                aligned_x1.append(i)

            if aa2 not in '-':
                aligned_x2.append(i)

        plots[0].x = np.array(aligned_x1)
        plots[1].x = np.array(aligned_x2)

        percent_id = sum(
                x[0] == x[1] and '-' not in x
                for x in zip(aligned_seq1, aligned_seq2)
        )
        percent_id /= max(len(p.seq) for p in plots)

        print(f"Scaffolds aligned with {100*percent_id:.2f}% identity.")

    if os.fork():
        sys.exit()

    for p in plots:
        plt.plot(p.x, p.y, label=p.label)

    plt.xlabel("aligned residue index" if args['--align'] else "residue index")
    plt.ylabel("relative deletions" if args['--normalize'] else "deletions")
    plt.xlim((0, max(p.x[-1] for x in plots)))
    plt.legend(loc='best')
    plt.show()

