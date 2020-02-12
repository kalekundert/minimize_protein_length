#!/usr/bin/env python3

import numpy as np
import pandas as pd

from pymol import cmd
from pymol.wizard import Wizard
from minp import MsaWorkspace, DeletionsWorkspace
from minp import load_weighted_msa, calc_deletion_scores, count_deletions

def paint_deletions(dels_workspace):
    """
DESCRIPTION

    Color the protein by how often residue is proposed to be deleted.

USAGE

    paint_deletions dels_workspace

ARGUMENTS

    dels_workspace = a path to a workspace for the final step of the pipeline, 
    which is to pick deletions.

NOTES

    This script assumes that the only structure loaded in pymol is the one 
    contained in the given workspace.
"""

    work_dels = DeletionsWorkspace.from_path(dels_workspace)
    msa = load_weighted_msa(work_dels.msa)
    dels = pd.read_hdf(work_dels.deletions_hdf5)
    n_dels = count_deletions(msa, dels)

    low = 0
    high = max(n_dels)

    sele = f'chain B and polymer.protein'
    cmd.alter(sele, 'b=n_dels[int(resi)-1]', space={**locals(), **globals()})
    cmd.spectrum('b', 'blue_white_red', sele, minimum=low, maximum=high)

    print(f'n:    {len(dels)}')
    print(f'low:  {low:.2f}')
    print(f'high: {high:.2f}')

def paint_deletion_scores(msa_workspace, low=0, high=50):
    """
DESCRIPTION

    Color each residue by its "deletion score", a measure of how commonly it is 
    deleted in homologous sequences.

USAGE

    paint_deletion_scores dels_workspace

ARGUMENTS

    dels_workspace = a path to a workspace for the final step of the pipeline, 
    which is to pick deletions.

    low = the percentile of the data to make the most blue {default: 0}.

    high = the percentile of the data to make the most red {default: 50}.

NOTES

    This script assumes that the only structure loaded in pymol is the one 
    contained in the given workspace.
"""
    work_blast, work_msa = MsaWorkspace.from_path(msa_workspace)
    msa = load_weighted_msa(work_msa)
    scores = calc_deletion_scores(msa)

    cutoff_low = np.percentile(scores, float(low))
    cutoff_high = np.percentile(scores, float(high))

    sele = f'chain B and polymer.protein'
    cmd.alter(sele, 'b=scores[int(resi)-1]', space={**locals(), **globals()})
    cmd.spectrum('b', 'blue_white_red', sele, minimum=cutoff_low, maximum=cutoff_high)

    print(f'low:  {cutoff_low:.2f} ({float(low):.2f}%)')
    print(f'high: {cutoff_high:.2f} ({float(high):.2f}%)')

def cycle_deletions(dels_workspace, cursor=0, low=0, high=50):
    """
DESCRIPTION

    Launch a wizard that will show each proposed deletion in turn.

USAGE

    cycle_delections dels_workspace [, cursor [, low [, high ]]]

ARGUMENTS

    dels_workspace = a path to a workspace for the final step of the pipeline, 
    which is to pick deletions.

    cursor = the deletion to start on {default: 0}

    low = the percentile of the data to make the most blue {default: 0}.

    high = the percentile of the data to make the most red {default: 50}.

NOTES

    This script assumes that the only structure loaded in pymol is the one 
    contained in the given workspace.
"""
    wizard = CycleDeletions(dels_workspace, cursor, low, high)
    cmd.set_wizard(wizard)

class CycleDeletions(Wizard):

    def __init__(self, dels_workspace, cursor=0, low=0, high=50):
        work_dels = DeletionsWorkspace.from_path(dels_workspace)
        msa = load_weighted_msa(work_dels.msa)
        scores = calc_deletion_scores(msa)

        self.dels = pd.read_hdf(work_dels.deletions_hdf5).sort_values(['del_start', 'del_end'])
        self.cursor = cursor
        self.low = float(low)
        self.high = float(high)
        self.sele = f'chain B and polymer.protein'

        cmd.alter(
                self.sele,
                'b=scores[int(resi)-1]',
                space={**locals(), **globals()},
        )
        self.redraw()

    def get_indices(self):
        i = int(self.dels.iloc[self.cursor]['del_start'] + 1)
        j = int(self.dels.iloc[self.cursor]['del_end'] + 1)
        return i, j

    def get_prompt(self):
        """
        Return text to be displayed in the top left corner of the view area.
        """
        i, j = self.get_indices()
        return [
                f"Showing deletion {self.cursor + 1}/{len(self.dels)}",
                f"resi {i}-{j}"
        ]

    def get_panel(self):
        """
        Return a list of lists describing the entries that should appear in the 
        right-hand panel of the GUI (below all the selections) when the wizard 
        is running.  The first few entries control common actions and basic 
        settings.  The remaining entries are buttons that the user can press to 
        see specific mutations.
        """

        # Each entry is described by a list with 3 elements.  The first element 
        # is a number that specifies the type of entry: 1 for text, 2 for a 
        # button, 3 for a menu.  The second element is the text that will be 
        # seen by the user.  The third element is an argument that means 
        # different things for each type of entry.  For buttons, it's a piece 
        # of code (given as a string) that should be executed when the button 
        # is pressed.  For menus, it's a "tag" that can be passed to get_menu()
        # to get a description of the menu in question.

        return [
            [1, "Wildtype vs Mutant Wizard", ''],
            [2, "Next", 'cmd.get_wizard().show_next()'],
            [2, "Previous", 'cmd.get_wizard().show_previous()'],
            [2, "Done", 'cmd.get_wizard().cleanup()'],
        ]

    def show_next(self):
        self.cursor = (self.cursor + 1) % len(self.dels)
        self.redraw()

    def show_previous(self):
        self.cursor = (self.cursor - 1) % len(self.dels)
        self.redraw()

    def redraw(self):
        cmd.refresh_wizard()

        i, j = self.get_indices()
        sele = f'({self.sele}) and resi {i}-{j}'
        print(sele)

        cmd.color('green', self.sele)
        cmd.spectrum('b', 'blue_white_red', sele, minimum=self.low, maximum=self.high)
        cmd.zoom(sele, buffer=10, animate=-1)

    def cleanup(self):
        cmd.set_wizard()

pymol.cmd.extend('paint_deletions', paint_deletions)
pymol.cmd.extend('paint_deletion_scores', paint_deletion_scores)
pymol.cmd.extend('cycle_deletions', cycle_deletions)





