#!/usr/bin/env python3

"""\
Download sequences with similarity to Cas9 from the NCBI databases

Usage:
    blast_cas9.py <workspace> [-d DATABASE] [-n HITS] [-f]

Arguments:
    <workspace>
        The directory where files relating to this run will be stored.  In 
        particular, this command will populate the workspace with the following 
        files:

        `blast_hits.xml`
            Results from the 
            and `blast_hits.fasta`
Options:
    -d --blast-database=NAME    [default: refseq_protein]
        The NCBI database to query.  The default is `refseq_protein`, which 
        excludes redundant sequences (important for Cas9, which has many exact 
        duplicates in the database).

    -n --blast-hits=NUM         [default: 5000]
        The number of BLAST hits to request.  Note that it's better to have too 
        many than too few, because the deletion-choosing algorithm will 
        naturally discount poor alignments.  But having more hits will also 
        make the alignment step more expensive.

    -f --force
        Clear any cached data.
"""

import docopt, toml
import sys, random, shlex
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from workspace import Workspace
from datetime import date
from pathlib import Path
from inform import display, fatal, plural

def int_or_none(x):
    return int(x) if x is not None else x

def read_meta(work):
    if work.blast_meta.exists():
        return toml.load(work.blast_meta)
    else:
        return dict(
                parameters=[],
                modification_dates={},
        )

def check_params(prev_params, curr_params):
    if prev_params and prev_params != curr_params:
        fatal("BLAST parameters differ from those used previously!")
        fatal()
        for key in params:
            prev = prev_params[key]
            curr = curr_params[key]
            if prev != curr:
                codicil(f"    {key} was {prev!r}, now {curr!r}")
        fatal()
        fatal("Use the -f flag to overwrite.  Aborting.")
        sys.exit(1)

args = docopt.docopt(__doc__)
work = Workspace(args['<workspace>'])
params = dict(
        query='WP_032461047.1',  # S. pyogenes Cas9
        database=args['--blast-database'],
        num_hits=int(args['--blast-hits']),
)
if args['--force']:
    work.rmdir()
    work.mkdir()

# Make sure the specified parameters are consistent with previous runs.
meta = read_meta(work)
check_params(meta['parameters'], params)
meta['parameters'] = params

# Query BLAST for Cas9 homologs.
if not work.blast_xml.exists():
    display("Performing BLAST query (30-45 min)...")
    response = NCBIWWW.qblast(
            'blastp', 
            params['database'],
            params['query'],
            hitlist_size=params['num_hits'],
            alignments=params['num_hits'],
    )

    work.blast_xml.write_text(response.read())
    meta['parameters'] = params
    meta['modification_dates'][work.blast_xml.name] = date.today()

# Download complete sequences for each BLAST hit.
if not work.blast_fasta.exists():
    display("Downloading hits from NCBI...")
    Entrez.email = 'kale_kundert@hms.harvard.edu'

    with work.blast_xml.open() as f:
        blast = NCBIXML.read(f)

    ids = [params['query']] + [
            x.hit_id.split('|')[1]
            for x in blast.alignments
    ]
    response = Entrez.efetch(db='protein', rettype='fasta', id=ids)

    work.blast_fasta.write_text(response.read())
    meta['modification_dates'][work.blast_fasta.name] = date.today()

seqs = list(SeqIO.parse(work.blast_fasta, 'fasta'))
display(f"Found {plural(seqs):# sequence/s}.")

meta['results'] = dict(num_hits=len(seqs))
with work.blast_meta.open('w') as f:
    toml.dump(meta, f)

