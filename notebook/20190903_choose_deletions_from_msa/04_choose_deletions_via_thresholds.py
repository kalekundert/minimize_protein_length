#!/usr/bin/env python3

"""\
Pick residues to delete by finding runs with above-average deletion scores.

Usage:
    04_choose_deletions_via_thresholds.py <msa_workspace>
"""

import numpy as np
import pandas as pd
from utils import *

def choose_deletions_via_thresholds(ref, msa):
    scores = calc_deletion_scores(ref, msa)
    dfs = []

    for threshold in scores:
        runs = find_runs(scores >= threshold)
        df = pd.DataFrame(runs, columns=['del_start', 'del_end'])
        dfs.append(df)

    dels = pd.concat(dfs)\
             .drop_duplicates()\
             .reset_index(drop=True)
    
    score_runs(dels, scores)
    dels.sort_values(by='score', ascending=False, inplace=True)

    # Require that each deletion score better than the average over the whole 
    # sequence.  This conveniently gets rid of the 0-N deletion.
    dels = dels[ dels['score'] > np.mean(scores) ]

    return dels.reset_index(drop=True)

def find_runs(a):
    # Taken from:
    # https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi

    # Make sure each run both starts and stops.
    bounded = np.hstack(([0], a, [0]))

    # This will give 1 at run starts and -1 at run ends.
    edges = np.diff(bounded)

    # The start and end indices refer to the positions before the runs starts 
    # and ends, respectively, which would be too small by 1.  But this error is 
    # inadvertently corrected by the fact that we added an element to the 
    # beginning of the input array, effectively making every index bigger by 1.  
    # So in the end, no further corrections are necessary.
    run_starts, = np.where(edges > 0)
    run_ends, = np.where(edges < 0)

    return zip(run_starts, run_ends)

def score_runs(dels, scores):
    dels['del_len'] = dels['del_end'] - dels['del_start']
    dels['score'] = dels.apply(
            lambda x: np.mean(scores[ int(x['del_start']) : int(x['del_end']) ]),
            axis=1,
    )


if __name__ == '__main__':

    # Parse the command-line arguments.

    import docopt
    args = docopt.docopt(__doc__)

    # Setup the workspaces

    work_blast, work_msa = MsaWorkspace.from_path(args['<msa_workspace>'])
    work_dels = DeletionsWorkspace(work_msa, 'del_threshold')

    # Calculate a list of candidate deletions.

    ref, msa = load_weighted_msa(work_msa)
    dels = choose_deletions_via_thresholds(ref, msa)

    # Record the results.

    print(f"Chose {len(dels)} deletions:\n")
    print(dels.describe())

    work_dels.write_deletions(ref, dels)
    work_dels.write_metadata()


def test_find_runs():
    examples = [
            ([], []),

            ([0], []),
            ([1], [(0, 1)]),

            ([0, 0], []),
            ([0, 1], [(1, 2)]),
            ([1, 0], [(0, 1)]),
            ([1, 1], [(0, 2)]),

            ([0, 0, 0, 0], []),
            ([1, 1, 1, 1], [(0, 4)]),
            ([1, 1, 0, 0], [(0, 2)]),
            ([0, 1, 1, 0], [(1, 3)]),
            ([0, 0, 1, 1], [(2, 4)]),
            ([1, 0, 0, 1], [(0, 1), (3, 4)]),
    ]

    for example, result in examples:
        assert list(find_runs(example)) == result

