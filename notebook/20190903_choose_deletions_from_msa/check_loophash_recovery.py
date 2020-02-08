#!/usr/bin/env python3

"""\
Test how well loophash recovers the width of known gaps in a structure.

Usage:
    03_check_loophash_recovery.py <msa_workspace> [-f]

Options:
    -f --force
        Overwrite cached data.
"""

import pandas as pd

from pyrosetta import init as init_rosetta, pose_from_file
from pyrosetta.rosetta.protocols import loophash
from pyrosetta.rosetta.utility import fixedsizearray1_double_6_t as Real6
from pyrosetta.rosetta.std import vector_unsigned_long
from utils import load_weighted_msa, MsaWorkspace, LoophashWorkspace

from importlib import import_module
choose_dels = import_module('04_choose_deletions_via_loophash')

def recover_known_gaps(db, aligned_pose, gap_size, radius):
    from more_itertools import stagger

    n = len(aligned_pose.ref_ungapped)
    records = []

    for i, j in stagger(range(n), (0, gap_size)):

        try:
            resi_n = aligned_pose.resis[i]
        except KeyError:
            continue

        try:
            resi_c = aligned_pose.resis[j]
        except KeyError:
            continue

        try:
            record = dict(
                    radius=radius,
                    i=i,
                    j=j,
                    resi_n=resi_n,
                    resi_c=resi_c,
                    smallest_hit=eval_segment(
                        db,
                        aligned_pose.pose,
                        resi_n,
                        resi_c,
                        radius,
                    )
            )
            records.append(record)

        except choose_dels.NoHitsFound:
            pass

    return pd.DataFrame(records)

def eval_segment(db, pose, resi_n, resi_c, radius):
    pose_bb = loophash.BackboneSegment()
    pose_bb.read_from_pose(pose, resi_n, resi_c - resi_n)

    # Find the 6D transformation between the ends of the gap.
    rt = Real6();
    success = loophash.get_rt_over_leap_without_foldtree_bs(
            pose, resi_n, resi_c, rt
    )
    if not success:
        raise No6dTransform

    # Find the shortest loops that can bridge the gap.
    hits = vector_unsigned_long()
    for k, lhm in sorted(db.lhms.items()):
        lhm.radial_lookup(radius, rt, hits)
        if len(hits) > 0:
            break
    else:
        raise NoHitsFound

    return k


if __name__ == '__main__':
    import docopt
    args = docopt.docopt(__doc__)

    work_blast, work_msa = MsaWorkspace.from_path(args['<msa_workspace>'])
    work_dels = LoophashWorkspace(work_msa)

    if not work_dels.recovery_hdf5.exists() or args['--force']:
        init_rosetta()

        msa = load_weighted_msa(work_msa)
        db = choose_dels.LoophashDatabases.load()
        aligned_pose = choose_dels.AlignedPose.from_file(
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

    df['len'] = df['j'] - df['i']
    df['recovered'] = (df['smallest_hit'] == df['len'])
    print(df.groupby(['radius'])['recovered'].mean())





