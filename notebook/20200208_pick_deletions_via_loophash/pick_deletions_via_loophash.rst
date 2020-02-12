***************************
Pick deletions via loophash
***************************

There are two issues with the algorithm detailed in 
:expt:`20200208_pick_deletions_via_threshold` that could both be addressed by 
considering the structure of the protein being shortened.  The first is that 
deletions spanning large gaps in 3D space are unlikely to yield a functional 
protein.  The second is that it may be too conservative to only delete 
contiguous runs of high-scoring residues (e.g. maybe more residues could be 
deleted to make a smaller gap in 3D space, or maybe it makes sense to also 
delete a few lower-scoring residues to make the ends of the deletion line up 
better).

An algorithm to determine whether a gap between two residues could easily be 
closed could address these issues.  A naive approach would be to simply measure 
the distance between the two residues.  This might work, but a more principled 
approach is to use loophash, and this approach has several advantages:

- The relative orientation of the takeoff and landing points is accounted for.

- The loophash hit provides a bound on how small the deletion can be.  No point 
  deleting 1 residue if you'll need to replace it with 2.

- The loophash hit could be used to replace the deleted region, making it more 
  likely that the adjacent structural elements can stay in place.  I'm not sure 
  if Xiao's assay can do this, though, and it might be better to keep the 
  native sequence anyways.


Algorithm
=========
- Enumerate all pairs of residues.

- Use loophash to find the smallest loop capable to connecting each pair, if 
  the residues between them were to be deleted.

- Require that the loop be smaller than a certain size (e.g. 6 residues) and 
  smaller that the residues it would replace.
   
- Let N = the number of residues between the pair minus the length of the 
  spanning loop.  This is the number of residues that will be deleted.

- Calculate the average deletion score for each possible deletion of length N 
  between the two residues.  Pick the window with the highest score to delete.

Methods
=======
::

   $ minp 04/loophash sp_cas9/best_msa
   All possible residue pairs:          935028           (1368 residues)
   → Residue missing from pose:         852165   −82863  (62 missing residues; expect −82863)
   → No spanning loop found:            222552  −629613
   → Gap smaller than spanning loop:    214304    −8248
   → Spanning loop >6 residues:          11289  −203015
   → Deleted >10% of the protein:         8909    −2380
   → Below-average deletion score:        4154    −4755
   → Duplicates:                          2428    −1726

            del_start      del_end    del_score    gap_start      gap_end      len_gap  len_spanning_loop  len_deletion
   count  2428.000000  2428.000000  2428.000000  2428.000000  2428.000000  2428.000000        2428.000000   2428.000000
   mean    826.897858   862.037891   108.572329   822.786540   863.477604    40.691063           5.551031     35.140033
   std     378.840959   390.065508    75.521138   378.860947   390.043688    35.063956           0.582220     34.995476
   min      10.000000    34.000000    45.965278     4.000000    33.500000     4.000000           3.000000      1.000000
   25%     530.000000   575.000000    59.039640   526.000000   576.875000    13.000000           5.000000      7.000000
   50%     927.000000   996.000000    82.065765   920.000000   996.000000    29.000000           6.000000     24.000000
   75%    1158.500000  1200.000000   128.550702  1158.000000  1201.500000    55.083333           6.000000     50.000000
   max    1341.000000  1353.000000   741.549708  1341.000000  1359.000000   142.000000           6.000000    136.000000

Results
=======

Radial lookup
-------------
The ``LoopHashMap::radial_lookup()`` is not really documented, so I had to read 
the code to see what it does.  I don't understand all of the code, but I'm 
pretty sure that ``radial_lookup()`` just expands the hash lookup to include 
adjacent bins.  So a radius of 1 will include one extra bin in each dimension, 
2 will include 2, etc.  I believe that a radius of 0 is equivalent to 
``lookup()``.

For this application, I want loophash to accurately represent how many residues 
are needed to span a gap between two residues.  It's not clear if a larger or 
smaller radius would be better for this, so I benchmarked the recovery of known 
gap sizes in spCas9.  Briefly, I made every possible deletion of 3–6 residues, 
used loophash with every radius from 0–4 to close the deletion, then calculated 
the fraction of the time that loophash chose a loop of the same length that was 
deleted::

   $ minp recovery sp_cas9/best_msa/

.. datatable:: loophash_recovery.csv

Loop size
---------
The most effective way to control the number of deletions proposed is to set a 
maximum length for the spanning loop.  Here this maximum is 6 residues.  This 
threshold eliminates ~200K potential deletions (probably ~50K of which would've 
passed the other filters).

I didn't really experiment with this parameter; I just view it as a way to 
control the number of results.

Deletions
---------
The loophash algorithm predicts about 4x more deletions than the threshold 
algorithm.  These are mostly focused on the same regions (e.g. those with high 
deletion scores):

.. figure:: deletions.svg

