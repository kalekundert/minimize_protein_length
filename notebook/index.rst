***********************
Minimize protein length
***********************

The goal of this project is to create an algorithm that predicts ways to 
shorten a protein without adversely affecting function.  Xiao has developed a 
high-throughput assay that can be used to test these deletions, so it is not 
necessary that very deletion be a good one.  Instead, we are just trying to 
create a pool of deletions enriched for those that will not affect function, to 
minimize the cost and effort of preparing/screening libraries.

The general approach we have envisioned is to use sequence alignments to 
identify residues that are potentially expendable.  There are 3 broad steps 
such an algorithm would need to complete:

1. Find homologous sequences.

2. Create a multiple sequence alignment (MSA).

3. Pick specific deletions to make.

.. todolist::

.. toctree::
   :maxdepth: 1
   :glob:
   :hidden:

   ideas_and_feedback.rst
   unexpected_observations.rst

   /20190808_find_homologs_via_blast/*
   /20190808_find_homologs_in_viral_genomes/*
   /20190815_find_homologs_via_structural_alignment/*
   /20200211_find_homologs_via_jackhmmer/*
   /20190903_compare_msa_algorithms/*
   /20200208_pick_deletions_via_threshold/*
   /20200208_pick_deletions_via_loophash/*
   /20200209_pick_deletions_for_cj_cas9/*
   /20200209_pick_deletions_for_sa_cas9/*
