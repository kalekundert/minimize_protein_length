#!/usr/bin/env python3

"""\
Download sequences with homology to the target.

Usage:
    minp 02_find_homologs <workspace> [-d DATABASE] [-n HITS] [-f]

Arguments:
    <workspace>
        The directory where files relating to this run will be stored.  In 
        particular, this command will populate the workspace with the following 
        files:

        `homologs.xml`
            Results from the BLAST query, in a parseable XML format.

        `homologs.fasta`
            Full sequences corresponding to each hit from the BLAST query.

        `metadata.toml`
            Information about the parameters used for the BLAST query.

Options:
    -d --database=NAME      [default: refseq_protein]
        The NCBI database to query.  The default is `refseq_protein`, which 
        excludes redundant sequences (important for Cas9, which has many exact 
        duplicates in the database).

    -n --num-hits=NUM       [default: 5000]
        The number of BLAST hits to request.  Note that it's better to have too 
        many than too few, because poor hits can be left out when building the 
        MSA, and some MSA algorithms may be better or worse than others at 
        aligning poor matches.

    -f --force
        Clear any cached data.
"""

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from inform import display, plural
from minp import HomologsWorkspace, report_elapsed_time
from io import StringIO

def main():

    # Parse the command-line arguments:

    import docopt
    args = docopt.docopt(__doc__)
    params = dict(
            algorithm='blast',
            database=args['--database'],
            num_hits=int(args['--num-hits']),
    )

    # Setup the workspace (e.g. where input and output files will go):

    work = HomologsWorkspace.from_params(
            args['<workspace>'],
            **params,
    )

    if args['--force']:
        work.rmdir()
        work.mkdir()

    work.update_params(**params)

    # Query BLAST for Cas9 homologs:

    if not work.output_xml.exists():
        display("Performing BLAST query (30-45 min)...")

        with report_elapsed_time():
            with work.touch(work.input_fasta) as p:
                SeqIO.write(work.shared.target_protein_record, p, 'fasta')

            with work.touch(work.output_xml) as p:
                response = NCBIWWW.qblast(
                        'blastp', 
                        work.database,
                        work.query,
                        hitlist_size=work.num_hits,
                        alignments=work.num_hits,
                )
                p.write_text(response.read())

    # Download complete sequences for each BLAST hit:

    if not work.output_fasta.exists():
        display("Downloading hits from NCBI...")
        Entrez.email = work.shared.user_email

        with report_elapsed_time():
            with work.output_xml.open() as f:
                blast = NCBIXML.read(f)

            ids = [
                    x.hit_id.split('|')[1]
                    for x in blast.alignments
            ]

            response = Entrez.efetch(db='protein', rettype='fasta', id=ids)
            with work.touch(work.output_fasta) as p:
                p.write_text(response.read(), p)

    # Report the results:

    seqs = list(SeqIO.parse(work.output_fasta, 'fasta'))
    display(f"Found {plural(seqs):# homolog/s}.")

    work.metadata['results'] = dict(num_hits=len(seqs))
    work.write_metadata()

