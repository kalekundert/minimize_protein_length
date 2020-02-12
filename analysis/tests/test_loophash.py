#!/usr/bin/env python3

import pytest

from pyrosetta import init as init_rosetta, pose_from_file
from minp.loophash import get_chain_info, map_pose_indices, pick_deletion_window

@pytest.fixture(scope='session')
def rosetta():
    init_rosetta()

def test_get_chain_info(rosetta):
    pose = pose_from_file('kysi.pdb')
    seq, resis = get_chain_info(pose, 'B')

    assert seq == 'KYSI'
    assert resis == [3,4,5,6]

@pytest.mark.parametrize(
        'path, chain, ref, expected', [(
            'kysi.pdb', 'B',
            'KYSI',
            {0:3, 1:4, 2:5, 3:6},
        ), (
            'kysi.pdb', 'B',
            'AKYSI',
            {1:3, 2:4, 3:5, 4:6},
        ), (
            'kysi.pdb', 'B',
            'KYASI',
            {0:3, 1:4, 3:5, 4:6},
        ), (
            'kysi.pdb', 'B',
            'KYSIA',
            {0:3, 1:4, 2:5, 3:6},
        ), (
            'kysi.pdb', 'B',
            'ASIA',
            {0:3, 1:4, 2:5, 3:6},
        )]
)
def test_map_pose_indices(path, chain, ref, expected, rosetta):
    pose = pose_from_file(path)
    mapping = map_pose_indices(ref, pose, chain)
    assert mapping == expected
    
@pytest.mark.parametrize(
        'scores, i, j, n, del_start, del_end, del_score', [(
            # Pick beginning (no sliding)
            [2, 4, 6, 8, 10],
             0, 1,            2,
             0, 2,            3,
        ), (
            # Pick end (no sliding)
            [2, 4, 6, 8, 10],
                      3,  4,  2,
                      3,  5,  9,
        ), (
            # Pick middle (no sliding)
            [2, 4, 6, 8, 10],
                1, 2,         2,
                1, 3,         5,
        ), (
            # Pick beginning (sliding)
            [10, 8, 6, 4, 2],
              0,    2,        2,
              0, 2,           9,
        ), (
            # Pick end (sliding)
            [2, 4, 6, 8, 10],
                   2,     4,  2,
                      3,  5,  9,
        )]
)
def test_pick_deletion_window(scores, i, j, n, del_start, del_end, del_score):
    expected = del_start, del_end, pytest.approx(del_score)
    assert pick_deletion_window(scores, i, j, n) == expected

