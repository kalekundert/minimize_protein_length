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
            Results from the BLAST query, in a parseable XML format.

        `blast_hits.fasta`
            Full sequences corresponding to each hit from the BLAST query.

        `metadata.toml`
            Information about the parameters used for the BLAST query.

Options:
    -d --blast-database=NAME    [default: refseq_protein]
        The NCBI database to query.  The default is `refseq_protein`, which 
        excludes redundant sequences (important for Cas9, which has many exact 
        duplicates in the database).

    -n --blast-hits=NUM         [default: 5000]
        The number of BLAST hits to request.  Note that it's better to have too 
        many than too few, because poor hits can be left out when building the 
        MSA, and some MSA algorithms may be better or worse than others at 
        aligning poor matches.

    -f --force
        Clear any cached data.
"""

import docopt
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from workspace import BlastWorkspace, report_elapsed_time
from inform import display, plural

if __name__ == '__main__':

    # Parse the command-line arguments:

    args = docopt.docopt(__doc__)
    params = dict(
            query='WP_032461047.1',  # S. pyogenes Cas9
            database=args['--blast-database'],
            num_hits=int(args['--blast-hits']),
    )

    # Setup the workspace (e.g. where input and output files will go):

    work = BlastWorkspace(args['<workspace>'])

    if args['--force']:
        work.rmdir()
        work.mkdir()

    work.update_params(**params)

    # Query BLAST for Cas9 homologs:

    if not work.output_xml.exists():
        display("Performing BLAST query (30-45 min)...")

        with report_elapsed_time():
            with work.touch(work.output_xml) as p:
                response = NCBIWWW.qblast(
                        'blastp', 
                        params['database'],
                        params['query'],
                        hitlist_size=params['num_hits'],
                        alignments=params['num_hits'],
                )
                p.write_text(response.read())

    # Download complete sequences for each BLAST hit:

    if not work.output_fasta.exists():
        display("Downloading hits from NCBI...")
        Entrez.email = 'kale_kundert@hms.harvard.edu'

        with report_elapsed_time():
            with work.output_xml.open() as f:
                blast = NCBIXML.read(f)

            ids = [
                    x.hit_id.split('|')[1]
                    for x in blast.alignments
            ]
            assert params['query'] in ids

            with work.touch(work.output_fasta) as p:
                response = Entrez.efetch(db='protein', rettype='fasta', id=ids)
                p.write_text(response.read())

    # Report the results:

    seqs = list(SeqIO.parse(work.output_fasta, 'fasta'))
    display(f"Found {plural(seqs):# sequence/s}.")

    work.metadata['results'] = dict(num_hits=len(seqs))
    work.write_metadata()

