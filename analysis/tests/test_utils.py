#!/usr/bin/env python3

import pytest
from minp.utils import *
from utils import msa_from_str, str_from_msa

@pytest.mark.parametrize(
        'seq, expected', [(
            'ABC',
            'ABC',
        ), (
            '.BC',
            'BC',
        ), (
            'A-C',
            'AC',
        ), (
            'AB.',
            'AB',
        )]
)
def test_remove_gaps(seq, expected):
    assert remove_gaps(seq) == expected

@pytest.mark.parametrize(
        'seq, expected', [(
            'AAA',
            {0:0, 1:1, 2:2},
        ), (
            '.AA',
            {1:0, 2:1},
        ), (
            'A-A',
            {0:0, 2:1},
        ), (
            'AA.',
            {0:0, 1:1},
        )]
)
def test_ungapped_indices(seq, expected):
    assert map_ungapped_indices(seq) == expected

@pytest.mark.parametrize(
        'msa_str, weights_expected', [(
            'M', [1],
        ), (
            'M/-', [1, 0],
        ), (
            'MN/M-', [1, 1/2],
        ), (
            'M-/MN', [1, 1],
        ), (
            'MN-P/M-QP', [1, 2/3],
        )]
)
def test_weight_alignments(msa_str, weights_expected):
    # Previously I was dividing by the length of the whole alignment, so I'd 
    # get wrong results if there were lots of gaps.  Make sure that's not still 
    # happening...
    
    msa = msa_from_str(msa_str)
    msa.ref = msa[0]

    weight_alignments(msa)
    weights_actual = [x.weight for x in msa]

    assert weights_actual == pytest.approx(weights_expected)

@pytest.mark.parametrize(
        'msa_str, scores_expected', [(
            'M', [0],
        ), (
            # no deletion
            'M/M', [0],
        ), (
            # no homology
            'M/-', [0],
        ), (
            'MN/M-', [0, 1/2],
        ), (
            'M-/MN', [0],
        ), (
            'MN-P/M-QP', [0, 2/3, 0],
        )]
)
def test_calc_deletion_scores(msa_str, scores_expected):
    msa = msa_from_str(msa_str)
    msa.ref = msa[0]

    weight_alignments(msa)
    scores_actual = calc_deletion_scores(msa)

    assert scores_actual == pytest.approx(scores_expected)
