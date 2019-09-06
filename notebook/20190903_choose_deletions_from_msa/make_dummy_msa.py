#!/usr/bin/env python3

from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalwCommandline
import random

cmd = ClustalwCommandline('clustalw2', infile='blast_refseq_100.fasta')
stdout, stderr = cmd()

msa = AlignIO.read('blast_refseq_100.aln', 'clustal')
print(msa)
