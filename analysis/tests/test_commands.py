#!/usr/bin/env python3

import pytest
from importlib import import_module
from utils import msa_from_str, str_from_msa

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

@pytest.mark.parametrize(
        'input, expected', [
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
)
def test_04_find_runs(input, expected):
    mod = import_module('minp.commands.04_choose_deletions_via_thresholds')
    assert list(mod.find_runs(input)) == expected

@pytest.mark.parametrize(
        'seqs_in, seqs_out', [(
            [],
            [],
        ), (
            [('1', 'MA')],
            [('1', 'MA')],
        ), (
            [('1', 'MA'), ('1', 'MA')],
            [('1', 'MA')],
        ), (
            [('1', 'MA'), ('2', 'MC')],
            [('1', 'MA'), ('2', 'MC')],
        )]
)
def test_03_remove_duplicates(seqs_in, seqs_out):
    mod = import_module('minp.commands.03_build_msa')

    def records_from_seqs(seqs):
        return [
                SeqRecord(
                    Seq(seq, IUPAC.protein),
                    id=id,
                )
                for id, seq in seqs
        ]

    def seqs_from_records(records):
        return [
                (x.id, str(x.seq))
                for x in records
        ]

    records_in = records_from_seqs(seqs_in)
    records_uniq = mod.remove_duplicates(records_in)
    seqs_uniq = seqs_from_records(records_uniq)
    
    assert seqs_uniq == seqs_out

@pytest.mark.parametrize(
        'ref_id,msa_in_str,msa_out_str', [(
            '1', "MNT", "MNT",
        ), (
            '1', "MNT/MNT", "MNT",
        ), (
            '1', "MNT/M-T", "MNT/M-T",
        ), (
            '1', "M-NT/MA-T", "MNT/M-T",
        )]
)
def test_03_prune_gaps(ref_id, msa_in_str, msa_out_str):
    mod = import_module('minp.commands.03_build_msa')

    msa_in = msa_from_str(msa_in_str)
    msa_pruned = mod.prune_gaps(msa_in, ref_id)
    msa_pruned_str = str_from_msa(msa_pruned)

    assert msa_pruned_str == msa_out_str
