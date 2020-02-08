*************************
Choose deletions from MSA
*************************

Once we have an MSA we're happy with, we need to use that MSA to propose 
reasonable deletions.  The purpose of this experiment is to play around with 
different ways of doing this.

To do
=====
- Strategy for terminal deletions.  The way I'm setting up my MSAs so far (i.e.  
  not counting terminal gaps) discourages terminal deletions.

Dummy MSA
=========
To have some data to work with, I created a quick MSA by:

- BLAST for 5000 protein sequences similar to Cas9
- Create an MSA for 100 of these sequences, randomly chosen.

This process can be reproduced (although the results will be pseudo-random) 
with the following scripts::

   ./blast_cas9.py
   ./make_dummy_msa.py

Reliability score threshold
===========================
The first part of this idea is to create a "reliability score" for each 
sequence in the alignment.  This score is meant to tell us, "if we see a 
deletion in this alignment, how seriously should we take that?"  A simple  
approach that I saw used in the literature for similar purposes is to calculate 
the percent identity of the alignment to the reference sequence (full-length 
Cas9).  For this calculation, I skipped:

- Terminal gaps: I don't want to penalize sequences that align really well with 
  only a certain domain, e.g. HNH.

- Gaps in the reference sequence: this would give more and more advantage to 
  sequences very closely related the reference as more and more distant 
  sequences are added to the alignment, creating more and more gaps in the MSA.  
  This works against part of the reason I want a reliability score, which is so 
  that the effect of adding bad sequences to the alignment will tend towards 
  zero (this will allow me to err on the side using more inclusive alignments).

The second part of the idea is to, for each ungapped position in the reference 
sequence, sum the reliability scores for each gap aligned to that position.  In 
other words, the more high quality alignments contain a gap at a certain 
position, the more confident we can be that we can delete that position.

We can see a plot of these sums:

.. figure:: blast_refseq_100_deletion_scores.svg

These scores seem to make sense when painted on the structure of Cas9.  Most of 
the high-scoring regions are on the surface and away from the nucleic acids.  
Some of the high scoring regions also match up well with places where a 
deletion wouldn't create a big gap in the 3D structure.  That said, I wouldn't 
read too much into this, because the underlying MSA I'm using is only based on 
100 sequences.

:download:`deletion_scores.pse`

The third and final part of this idea is that we can set a threshold (i.e. a 
vertical line on the above plot) and delete any contiguous stretch of residues 
that score higher than that.  Then, we can slide this threshold all the way 
from the highest score to the lowest to create a large set of plausible 
deletions.  Deletions that score worse than the average score over the whole 
protein are discarded.  This is not a stringent threshold, but it prevents 
proposing that we should delete the whole protein or small groups of residues 
based purely on bad alignments.

The above scheme is implemented by the following script::

   ./choose_deletions_via_thresholds.py

Small gaps
==========
One issue I can foresee with the above approach is that deletions spanning 
large gaps in 3D space are unlikely to yield a functional protein.  To account 
for this, I came up with the following algorithm:

- Enumerate all possible deletions (i.e. all combinations of residues)

- Discard those that would create ends more than 15Å apart.  I chose 15Å by 
  looking at some distances between adjacent secondary structures in pymol, but 
  it's not really anything quantitative.

- Score the remaining deletions as above, and discard any that would score 
  worse than the average over the whole structure.

This scheme is implemented by the following script::

   ./choose_deletions_via_gaps.py

This approach produces more deletions, about ~10,000 compared to ~500 for the 
threshold-based algorithm.  I'm sure that most of this comes from just more 
finely subdividing regions that generally score highly.

Loophash
========
While the above strategy simply uses an arbitrary cutoff to how far is too far, 
a more principled approach is to use loophash (perhaps limited to 4 residues) 
to attempt to close the gap.  This provides two further advantages:

- The loophash hit can be used to replace the deleted region, making it more 
  likely that the adjacent structural elements can stay in place.

- The loophash hit also provides a bound on how small the deletion can be.  No 
  point deleting 1 residue if you're going to replace it with 2.

I'm not familiar enough with loophash to know how difficult this will be, or 
how much performance will be an issue.




BLAST
=====
The BLAST result are fine, I checked that I get the same top sequences via the 
web interface.  But somehow I'm getting the wrong sequence for the query 
sequence...  I think it has to do with putting the query sequence in the MSA 
twice.  Fixing that seems to fix the really bad alignments.

MSA
===
I like the idea of evaluating the MSA algorithms and determining an appropriate 
number of input sequences based on the number of aligned sequences with >30% 
identity to the query.  If too few sequences are given, this value will not be 
maximized.  If too many sequences are given, the alignments seem to start 
failing, and this value goes down again.  If the alignment algorithm is too 
slow to practically maximize this value, then it's not a viable candidate.

.. datatable:: msa_times.xlsx

   Times are upper limits; I was often running multiple simulations at once 
   with the CPU maxed out.

A slightly improved version of the above scheme is to plot percent identity 
histograms for each MSA.  This shows not only which MSAs have to most sequences 
with >30% identity, but also which algorithms have higher identity in general.  
More instance, this reveals that Clustalω performs significantly better than 
MAFFT with 900 sequences, even though both align ≈870 sequences with >30% 
identity.

.. image:: mafft_vs_clustalo.svg

Using the approach, it currently seems like Clustalω with 900 sequences is the 
best alignment.

I'd be interested in actually using an optimization algorithm to figure out 
exactly how many sequences to use, but that's probably overkill for now; I know 
900 is pretty good.  Maybe this would be something for if this works and we 
want to distribute it.

Scoring
=======
Scaling the sequence weights to ignore alignments with <30% identity has a big 
effect on the results, even though it doesn't discard that many sequences.  I 
think the bigger effect is that it effectively up-weights the better aligned 
sequences.  I don't need to do that; those gaps score the highest anyways.  

Really what I'm concerned about is including complete garbage in the alignment.  
I mean, the garbage is there; I'm not being choosy in the BLAST or MSA steps.  
I just don't want to read too much into it.  Setting this 30% identity 
threshold does really get rid of it, but 

Loophash
========
- I think radius basically just checks adjacent 6D bins for hits.

- It looks like radius=1 gives the most faithful recovery of know gap length. 


