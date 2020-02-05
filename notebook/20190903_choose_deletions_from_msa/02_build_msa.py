#!/usr/bin/env python3

"""\
Create a MSA for sequences identified by BLAST and downloaded from NCBI.

Usage:
    build_msa <workspace> [-a ALGORITHM] [-l LIMIT] [-f]

Arguments:
    <workspace>
        The directory containing the sequences to align.

Options:
    -a --algorithm=NAME     [default: mafft]
        The name of the multiple sequence alignment algorithm to use.  The 
        following names are supported:

            mafft       muscle      clustalo    dialign2
            msaprobs    probcons    tcoffee

    -l --limit=NUM 
        Only include the given number of sequences in the alignment.  The query 
        sequence will always be included, but the others will be picked 
        randomly.

    -f --force
        Clear existing results (if any) and recalculate the MSA.
"""

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.Align import Applications, MultipleSeqAlignment
from utils import BlastWorkspace, MsaWorkspace, report_elapsed_time
from inform import display

def mafft_wrapper(work_msa):
    app = Applications.MafftCommandline(
            input=work_msa.input_fasta,
            clustalout=True,
    )
    stdout, stderr = app()
    work_msa.output_aln.write_text(stdout)

def muscle_wrapper(work_msa):
    app = Applications.MuscleCommandline(
            input=work_msa.input_fasta,
            out=work_msa.output_aln,
            clw=True,
    )
    app()

def clustalo_wrapper(work_msa):
    app = Applications.ClustalOmegaCommandline(
            infile=work_msa.input_fasta,
            outfile=work_msa.output_aln,
            outfmt='clustal',
            verbose=True,
            auto=True,
    )
    app()

def dialign2_wrapper(work_msa):
    raise NotImplementedError("I can't figure out how to get dialign to output the MSA it calculates...")
    app = Applications.DialignCommandline(
            'dialign',
            input=work_msa.input_fasta,
            fn=work_msa.output_aln.stem,
    )
    app()

def msaprobs_wrapper(work_msa):
    app = Applications.MSAProbsCommandline(
            infile=work_msa.input_fasta,
            outfile=work_msa.output_aln,
            clustalw=True,
    )
    app()

def probcons_wrapper(work_msa):
    app = Applications.ProbconsCommandline(
            input=work_msa.input_fasta,
            clustalw=True,
    )
    stdout, stderr = app()
    work_msa.output_aln.write_text(stdout)

def tcoffee_wrapper(work_msa):
    app = Applications.TCoffeeCommandline(
            infile=work_msa.input_fasta,
            outfile=work_msa.output_aln,
            output='clustalw',
    )
    app()

apps = {
        'mafft':    mafft_wrapper,
        'muscle':   muscle_wrapper,
        'clustalo': clustalo_wrapper,
        'dialign2': dialign2_wrapper,
        'msaprobs': msaprobs_wrapper,
        'probcons': probcons_wrapper,
        'tcoffee':  tcoffee_wrapper
}

def prune_gaps(msa, ref_id):
    """
    Remove from the given MSA:

    - Any positions that have gaps in the reference sequence.
    - Any sequences which don't have any gaps relative to the reference.

    This is meant to make it easier visually inspect the MSA and find gaps 
    that might represent residues that can be safely deleted.
    """

    # Find the reference sequence.
    for record in msa:
        if record.id == ref_id:
            ref = record
            break
    else:
        raise ValueError("reference sequence '{ref_id}' not found in MSA")

    # Find the indices that correspond to non-gaps in the reference sequence.
    non_gaps = []
    for i, nuc in enumerate(ref.seq):
        if nuc != '-':
            non_gaps.append(i)

    msa_pruned = MultipleSeqAlignment([])

    # Construct a new MSA with the pruned sequences.
    for record in msa:
        record.seq = Seq(
                ''.join(record.seq[i] for i in non_gaps),
                record.seq.alphabet,
        )

        if record.id == ref_id or '-' in record.seq:
            msa_pruned.append(record)

    return msa_pruned


if __name__ == '__main__':

    # Parse the command-line arugments:

    import docopt
    args = docopt.docopt(__doc__)
    params = dict(
        algorithm=args['--algorithm'],
        limit=int(args['--limit'] or 0),
    )

    # Setup the workspace (e.g. where input and output files will go):

    work_blast = BlastWorkspace(args['<workspace>'])
    work_msa = MsaWorkspace.from_params(work_blast, **params)

    if args['--force']:
        work_msa.rmdir()
        work_msa.mkdir()

    work_msa.update_params(**params)

    # Create an input file containing a subset of the BLAST hits (if requested):

    if not work_msa.input_fasta.exists():
        seqs = list(SeqIO.parse(work_blast.output_fasta, 'fasta'))
        if params['limit']:
            seqs = seqs[:params['limit']]

        with work_msa.touch(work_msa.input_fasta) as p:
            SeqIO.write(seqs, p, 'fasta')

    # Create the MSA:

    if not work_msa.output_aln.exists():
        display(f"Running '{params['algorithm']}'...")

        with report_elapsed_time():
            with work_msa.touch(work_msa.output_aln):
                apps[params['algorithm']](work_msa)

    # Create a simplified MSA that's easier to visually inspect:

    msa = AlignIO.read(work_msa.output_aln, 'clustal')
    msa_pruned = prune_gaps(msa, work_blast.query)

    if not work_msa.pruned_aln.exists():
        with work_msa.touch(work_msa.pruned_aln) as p:
            with p.open('w') as f:
                AlignIO.write(msa_pruned, f, 'clustal')

    display(msa)
    work_msa.write_metadata()

