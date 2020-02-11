#!/usr/bin/env python3

"""\
Setup a workspace containing the input files necessary to predict the best ways 
to shorten a target protein.

Usage:
    minp 01_init <workspace> <fasta> <pdb> <pdb_chain> <loophash_db> [options]

Arguments:
    <workspace>
        The path to a directory that will contain the input files.

    <fasta>
        The path to a *.fasta file containing the DNA sequence of the protein 
        to be shortened.  This is used to find homologous sequences, and to 
        report the specific sequences to delete.

    <pdb>
        The path to a structure of the protein to be shortened, in either the 
        PDB or mmCIF format.  This is used when using loophash to pick 
        deletions.

    <pdb_chain>
        The chain in the PDB file containing the protein of interest.

    <loophash_db>
        The path to a directory containing databases of loop fragments for use 
        with loophash.

Options:
    -f --force
        If the given workspace already exists, overwrite it.
        
    -e --user-email <ADDRESS>
        Your email address.  This is used to prevent throttling when submitting 
        requests to web servers such as BLAST.
"""

import os.path
from minp import SharedWorkspace
from shutil import copyfile
from pathlib import Path
from inform import display, error, fatal

def main():

    import docopt
    args = docopt.docopt(__doc__)

    # Create the workspace directory:

    work = SharedWorkspace(args['<workspace>'])

    if args['--force']:
        work.rmdir()

    if work.exists():
        error(f"Workspace '{args['<workspace>']}' already exists.")
        fatal("Use the -f flag to overwrite.  Aborting.")

    work.mkdir()

    # Fill in the workspace:

    copyfile(args['<fasta>'], work.target_fasta)
    copyfile(args['<pdb>'], work.get_target_pdb(args['<pdb>']))

    work.loophash_db.symlink_to(
            os.path.relpath(
                Path(args['<loophash_db>']).resolve(),
                work.root,
            ),
            target_is_directory=True,
    )

    work.params['target'] = {
            'pdb_chain': args['<pdb_chain>'],
    }
    work.params['user'] = {
            'email': args['--user-email'],
    }
    work.write_params()

    display("Workspace successfully initialized.")




