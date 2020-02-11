#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped

def msa_from_str(seqs_str):
    records = []

    for i, line in enumerate(seqs_str.split('/'), 1):
        seq = Seq(line, Gapped(IUPAC.protein, '-'))
        record = SeqRecord(seq, id=str(i))
        records.append(record)

    return MultipleSeqAlignment(records)

def str_from_msa(msa):
    return '/'.join(
            str(x.seq)
            for x in msa
    )
