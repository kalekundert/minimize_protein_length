#!/usr/bin/env python3

"""\
Test how well loophash recovers the width of known gaps in a structure.

Usage:
    minp check_loophash_recovery <msa_workspace> [-f]

Options:
    -f --force
        Overwrite cached data.
"""

import pandas as pd

from pyrosetta import init as init_rosetta
from minp import load_weighted_msa, MsaWorkspace, LoophashWorkspace
from minp.loophash import LoophashDatabases, AlignedPose

def main():

    # Setup the workspaces:

    import docopt
    args = docopt.docopt(__doc__)

    work_blast, work_msa = MsaWorkspace.from_path(args['<msa_workspace>'])
    work_dels = LoophashWorkspace(work_msa)

    # Calculate the size of the loops returned by loophash for known gaps:

    if not work_dels.recovery_hdf5.exists() or args['--force']:
        init_rosetta()

        msa = load_weighted_msa(work_msa)
        db = LoophashDatabases.load()
        aligned_pose = AlignedPose.from_file(
                msa,
                work_blast.input_pdb,
                work_blast.input_pdb_chain,
        )

        df = pd.concat([
            recover_known_gaps(db, aligned_pose, gap, radius)
            for gap in [3,4,5,6]
            for radius in [0,1,2,3,4]
        ])

        with work_dels.touch(work_dels.recovery_hdf5) as p:
            df.to_hdf(p, 'df')

    else:
        df = pd.read_hdf(work_dels.recovery_hdf5, 'df')

    # Display the results:

    df['len'] = df['j'] - df['i']
    df['recovered'] = (df['smallest_hit'] == df['len'])
    print(df.groupby(['radius'])['recovered'].mean())

def recover_known_gaps(db, aligned_pose, gap_size, radius):
    from more_itertools import stagger

    n = len(aligned_pose.ref_ungapped)
    records = []

    for i, j in stagger(range(n), (0, gap_size)):
        if not aligned_pose.has_resis(i, j):
            continue

        smallest_hit = find_smallest_spanning_loop(
                    db,
                    aligned_pose.pose,
                    aligned_pose.resis[i],
                    aligned_pose.resis[j],
                    radius,
        )
        record = dict(
                radius=radius,
                i=i,
                j=j,
                resi_n=aligned_pose.resis[i],
                resi_c=aligned_pose.resis[j],
                smallest_hit=smallest_hit,
        )
        records.append(record)

    return pd.DataFrame(records)
