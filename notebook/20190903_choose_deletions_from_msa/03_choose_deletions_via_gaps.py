#!/usr/bin/env python3

from choose_deletions_via_thresholds import *

import prody
from Bio.Data.IUPACData import protein_letters_3to1

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools as it

# Another idea:
#
# - Iterate through structure to find all possible segments with nearby takeoff 
#   and landings points.
#   
#   - Could use loophash to decide if endpoints are "close enough". 
#
#     Advantages:
#
#       - The identified loop could be used to replace the deleted residues.  
#         In turn, this use of canonical loops might help folding kinetics.
#
#       - The identified loop would also provide a natural limit on the size of 
#         deletions to consider, e.g. no reason to consider a 2-residue 
#         deletion if it would be replaced by a 4-residue turn.
#   
#       - Shouldn't be too expensive.
#
#     Disadvantages:
#
#       - Some good deletions might need to traverse fairly long gaps, e.g.  
#         with a Gly-Ser linker.  
#
#       - Not especially tunable.  Parameters really depend on what's in the 
#         loophash database.
#
#   - Could use a distance cutoff to determine what's close.
#
# - Score each potential deletion based on the multiple sequence alignment:
#
#   - Deletion score:
#
#       - How often is the position deleted in the MSA?
#
#       - For each sequence in MSA: weight = percent identity to reference.  
#         Should probably ignore terminal deletions.
#
#       - For each position in segment: score = summed weights of all positions 
#         that are deletions.
#
#   - Conservation score:
#
#       - How well is the position conserved in the MSA?
#
#       - Weight: same as deletion score.
#
#       - For each position in segment: score = BLOSUM/PAM (scaled to 1) of 
#         reference vs MSA sequence times weight of said sequence.
#
#   - For each deletion being considered: subtract the conservation score from 
#     the deletion score.
#     
#   - Rank the deletions by their combined scores.
#     

def load_ref_struct(pdb_path, sele):
    atoms = prody.parseCIF(pdb_path).select(sele)

    if atoms is None:
        raise ValueError(f"no atoms were selected by '{sele}'")

    return atoms

def choose_deletions_via_gaps(ref_seq, ref_struct, msa):
    from scipy.spatial.distance import euclidean

    ref_record, ref_seq = ref_seq, remove_gaps(ref_seq)
    ref_struct = ref_struct.getHierView()

    # Make sure the structure seems to match up with the sequences that the 
    # scores are based on.  Print out any mismatches so the user is aware, but 
    # don't require the sequences to be identical.
    
    if ref_struct.numChains() != 1:
        raise ValueError("can only design deletions for single-chain proteins.")

    mutations = []
    for res in ref_struct.iterResidues():
        aa_seq = ref_seq[res.getResnum() - 1]
        aa_struct = protein_letters_3to1[res.getResname().title()]

        if aa_seq != aa_struct:
            mutations.append(f'{aa_seq}{res.getResnum()}{aa_struct}')

    if mutations:
        from textwrap import fill, indent
        n = max(len(x) for x in mutations)
        print(f"Found {len(mutations)} differences ({100 * len(mutations) / ref_struct.numResidues():.1f}%) between the structure and the sequence:")
        print(indent(fill(" ".join(f'{x:{n}s}' for x in mutations)), "    "))
        print()

    # Iterate through every segment that could be deleted.  For those where the 
    # resulting chain ends would be nearby, use the MSA to predict how well- 
    # tolerated that deletion is likely to be.

    hits = []
    hit_threshold_A = 15  # Chosen by on looking at things in pymol.
    n_res = ref_struct.numResidues()
    n_combos = n_res * (n_res - 1) // 2
    progress = 1

    for res_n, res_c in it.combinations(ref_struct.iterResidues(), 2):
        sys.stdout.write(f'\r[{progress}/{n_combos}]')
        progress += 1

        resi_n = res_n.getResnum()
        resi_c = res_c.getResnum()
        assert resi_n < resi_c

        # How big a gap would this deletion create?
        xyz_n = res_n['C'].getCoords()
        xyz_c = res_c['N'].getCoords()
        gap = euclidean(xyz_n, xyz_c)

        if gap < hit_threshold_A:
            hits += [{
                'del_start': resi_n - 1,  # slice indexing
                'del_end': resi_c,
                'gap': gap,
            }]

    dels = pd.DataFrame(hits)

    scores = calc_deletion_scores(ref_record, msa)
    score_deletions(dels, scores)

    dels.sort_values(by='score', ascending=False, inplace=True)

    # Require that each deletion score better than average.  Not meant to be a 
    # stringent threshold, just meant to get rid of things that really don't 
    # make sense.
    dels = dels[ dels['score'] > np.mean(scores) ]

    # Show the user what happened.

    n = len(str(n_combos))
    print(f"\rFound {n_combos:{n}} possible deletions (from a model with {n_res} residues)")
    print(f"Found {len(hits):{n}} deletions with endpoints within {hit_threshold_A}Å")
    print(f"Kept  {len(dels):{n}} deletions with better-than-random MSA scores")

    return dels.reset_index(drop=True)

def choose_deletions_via_loophash(ref_seq, ref_struct, msa):
    pass


if __name__ == '__main__':
    msa_path = Path('blast_refseq_100.aln')
    pdb_path = Path('4un3.cif')

    ref_seq, msa = load_msa(msa_path, 'clustal', 'WP_032461047.1')
    ref_struct = load_ref_struct(pdb_path, 'protein and chain B')

    weight_alignments(ref_seq, msa)

    # Output a list of candidate deletions.
    dels = choose_deletions_via_gaps(ref_seq, ref_struct, msa)
    output_deletions(ref_seq, dels, f'{msa_path.stem}_gap_dels.fasta')
