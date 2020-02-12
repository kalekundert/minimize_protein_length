#!/usr/bin/env python3

"""\
Pick residues to delete by finding gaps that can be closed by loophash.

Usage:
    minp 04_choose_deletions_via_loophash <msa_workspace> [-f]

Options:
    -f --force
        Discard and recalculate any cached data.
"""

import numpy as np
import pandas as pd

from minp import MsaWorkspace, LoophashWorkspace
from minp import load_weighted_msa, calc_deletion_scores
from minp.loophash import load_spannable_gaps, choose_gaps_to_delete, GapFilters

def main():

    # Setup the workspaces:

    import docopt
    args = docopt.docopt(__doc__)

    work_blast, work_msa = MsaWorkspace.from_path(args['<msa_workspace>'])
    work_dels = LoophashWorkspace(work_msa)

    if args['--force']:
        work_dels.rmdir()
        work_dels.mkdir()

    # Choose which deletions to make:

    if work_dels.deletions_hdf5.exists():
        dels = pd.read_hdf(work_dels.deletions_hdf5)
        filters = GapFilters.from_toml(work_dels.filters_toml)
    else:
        msa = load_weighted_msa(work_msa)
        scores = calc_deletion_scores(msa)
        gaps, filters = load_spannable_gaps(work_dels, msa)
        dels = choose_gaps_to_delete(scores, gaps, filters)

    # Analyze the results:

    print(filters)
    print(dels.describe()); print()

    work_dels.write_deletions(dels, filters)
    work_dels.write_metadata()

