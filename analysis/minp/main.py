#!/usr/bin/env python3

"""\
Predict ways to shorten the length of a protein without adversely impacting 
function, using information from multiple sequence alignments.

Usage:
    minp <command> [<args>...]
    minp --version
    minp --help

Arguments:
    <command>
        The name of the command you want to run.  You only need to specify 
        enough of the name to be unique.  Broadly speaking, there are two 
        categories of scripts.  The first are part of the main design pipeline.
        These are prefixed with numbers so that you know the order to run them 
        in.  The second are helper scripts and are not prefixed.

{command_table}

    <args>...
        The necessary arguments depend on the command being run.  For more 
        information, pass the '--help' flag to the command you want to run.

Options:
    -v, --version
        Display the version of `minp` that's installed.

    -h, --help
        Display this help message.
"""

import sys, re, shlex
import docopt
from inform import fatal, Error
from itertools import zip_longest
from pathlib import Path
from . import __version__

def cmds_from_dir(dir):
    return [
            (path.stem, path)
            for path in dir.glob('*.py')
            if not path.name.startswith('__')
    ]

def make_command_table(cmds):
    """
    Return a nicely formatted table of all the `minp` commands installed on 
    this system to incorporate into the help text.  The table will have two 
    columns.  The first will list the commands that comprise the main pipeline 
    and the second will list all the other miscellaneous helper functions.
    """

    # Split every command installed on the system into two categories: those 
    # that are part of the main pipeline and those that are just utilities or 
    # helpers.  Pipeline scripts start with numbers, helper scripts don't.

    pipeline_cmds = []
    helper_cmds = []

    for cmd, main in sorted(cmds):
        if re.match('\d+_', cmd):
            pipeline_cmds.append(cmd)
        else:
            helper_cmds.append(cmd)

    # Make the table.

    rows = []
    columns = zip_longest(pipeline_cmds, helper_cmds, fillvalue='')
    longest_pipeline_cmd = max(len(x) for x in pipeline_cmds)

    for cols in columns:
        row = '        {0[0]:{1}}   {0[1]}'.format(cols, longest_pipeline_cmd)
        rows.append(row)

    return '\n'.join(rows)

def find_matching_command(given_cmd, known_cmds):
    matching_cmds = find_matching_commands(given_cmd, known_cmds)
    pseudo_cmd = shlex.join([Path(sys.argv[0]).name] + sys.argv[1:])


    if len(matching_cmds) == 0:
        guessed_cmd = did_you_mean(given_cmd, known_cmds)
        fatal(f"""\
Unknown command '{given_cmd}'.  Did you mean:

    $ {pseudo_cmd.replace(given_cmd, guessed_cmd)}

""")

    elif len(matching_cmds) > 1:
        err = f"Command '{given_cmd}' is ambiguous.  Did you mean:\n\n"
        for cmd, main in matching_cmds:
            err += f"    $ {pseudo_cmd.replace(given_cmd, cmd)}\n"
        err += '\n'
        fatal(err)

    else:
        return matching_cmds[0]

def find_matching_commands(given_cmd, known_cmds):
    from fnmatch import translate

    matching_cmds = []
    cmd_pattern = given_cmd.replace('/', '.*')

    for cmd, path in known_cmds:
        if re.search(cmd_pattern, cmd):
            matching_cmds.append((cmd, path))

    return matching_cmds

def update_sys_argv(given_cmd, full_cmd):
    """
    Update sys.argv to contain the full subcommand name, so that the argument 
    parser for the subcommand won't get confused.
    """
    i = sys.argv[1:].index(given_cmd)
    sys.argv[i+1] = full_cmd

def did_you_mean(unknown_cmd, known_cmds):
    """
    Return the command with the name most similar to what the user typed.  This 
    is used to suggest a correct command when the user types an illegal 
    command.
    """
    from difflib import SequenceMatcher
    similarity = lambda x: SequenceMatcher(None, x[0], unknown_cmd).ratio()
    did_you_mean = sorted(known_cmds, key=similarity, reverse=True)
    return did_you_mean[0][0]


def main():
    try:
        # Find all available commands:
        cmds = cmds_from_dir(Path(__file__).parent / 'commands')

        # Read the command the user typed on the command line:
        command_table = make_command_table(cmds)
        args = docopt.docopt(
                __doc__.format(**locals()),
                version=__version__,
                options_first=True,
        )
        given_cmd = args['<command>']

        # Execute the specified command:
        cmd, path = find_matching_command(given_cmd, cmds)
        update_sys_argv(given_cmd, cmd)

        from runpy import run_path
        module = run_path(str(path))
        module['main']()

    except Error:
        e.report()
