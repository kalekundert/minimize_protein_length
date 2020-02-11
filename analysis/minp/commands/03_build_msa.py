#!/usr/bin/env python3

"""\
Create a MSA for sequences identified by BLAST and downloaded from NCBI.

Usage:
    minp 03_build_msa <homologs_workspace> [-a ALGORITHM] [-l LIMIT] [-f]

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
from minp import MsaWorkspace, report_elapsed_time
from inform import display, plural

def main():

    # Parse the command-line arguments:

    import docopt
    args = docopt.docopt(__doc__)
    params = dict(
        algorithm=args['--algorithm'],
        limit=int(args['--limit'] or 0),
    )

    # Setup the workspace (e.g. where input and output files will go):

    work_msa = MsaWorkspace.from_params(
            args['<homologs_workspace>'],
            **params,
    )

    if args['--force']:
        work_msa.rmdir()
        work_msa.mkdir()

    work_msa.update_params(**params)

    # Create an input file containing the target sequence and (possibly) a 
    # subset of the homologs:

    if not work_msa.input_fasta.exists():
        seqs = remove_duplicates([
                work_msa.shared.target_protein_record,
                *work_msa.homologs.output_records,
        ])
        if work_msa.limit:
            seqs = seqs[:work_msa.limit]

        with work_msa.touch(work_msa.input_fasta) as p:
            SeqIO.write(seqs, p, 'fasta')

    # Create the MSA:

    if not work_msa.output_aln.exists():
        display(f"Running '{work_msa.algorithm}'...")

        with report_elapsed_time():
            with work_msa.touch(work_msa.output_aln):
                apps[work_msa.algorithm](work_msa)

    # Create a simplified MSA that's easier to visually inspect:

    msa = AlignIO.read(work_msa.output_aln, 'clustal')
    msa_pruned = prune_gaps(msa, work_msa.shared.target_id)

    if not work_msa.pruned_aln.exists():
        with work_msa.touch(work_msa.pruned_aln) as p:
            with p.open('w') as f:
                AlignIO.write(msa_pruned, f, 'clustal')

    display(msa)
    work_msa.write_metadata()


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

def remove_duplicates(seqs):
    already_seen = set()
    num_duplicates = 0
    unique_seqs = []

    for seq in seqs:
        if seq.id in already_seen:
            num_duplicates += 1
            continue

        already_seen.add(seq.id)
        unique_seqs.append(seq)

    if num_duplicates:
        display(f"Removed {plural(num_duplicates):# duplicate sequence/s}.")

    return unique_seqs

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
        raise ValueError(f"reference sequence '{ref_id}' not found in MSA")

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

