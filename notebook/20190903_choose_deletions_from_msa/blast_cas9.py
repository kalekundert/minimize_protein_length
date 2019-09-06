#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from pathlib import Path
import random

blast_xml = Path('blast_refseq.xml')
blast_fa = Path('blast_refseq.fasta')
cas9_id = 'WP_032461047.1'

if not blast_xml.exists():
    print("Performing BLAST query...")
    response = NCBIWWW.qblast(
            'blastp', 'refseq_protein', cas9_id,
            alignments=5000,
    )
    with blast_xml.open('w') as f:
        f.write(response.read())

with blast_xml.open() as f:
    blast = NCBIXML.read(f)

if not blast_fa.exists():
    print("Downloading hits from NCBI...")
    Entrez.email = 'kale_kundert@hms.harvard.edu'

    ids = [cas9_id] + [x.hit_id.split('|')[1] for x in blast.alignments]
    response = Entrez.efetch(db='protein', rettype='fasta', id=ids)

    with blast_fa.open('w') as f:
        f.write(response.read())

seqs = list(SeqIO.parse(blast_fa, 'fasta'))
seqs_100 = [seqs[0]] + random.sample(seqs[1:], 100)
SeqIO.write(seqs_100, f'{blast_fa.stem}_100{blast_fa.suffix}', 'fasta')

