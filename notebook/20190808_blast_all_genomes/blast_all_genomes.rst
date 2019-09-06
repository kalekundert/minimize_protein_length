*****************
BLAST all genomes
*****************

Methods
=======

BLAST
-----
Cas9 accession number::

   WP_032461047.1

I've found that the best way to run BLAST is through the biopython interface.  
Obviously using the web interface isn't very reproducible, although it's good 
for playing around with ideas.  BLAST has a command-line interface, but 
biopython is needed to parse the output anyways.  It might make sense to use 
the BLAST command-line interface to run local queries on O2, which has all of 
the NCBI BLAST databases, but they seem to be fairly out of date.  

I chose to use ``refseq_protein`` because ``nr`` has:

- A lot of nearly identical sequences from many species (presumably the ones 
  that have been sequenced multiple times).
- A lot of plasmids.

These sequences don't add any information about which residues we can delete, 
so if I were to use ``nr``, I'd want to filter them out before doing any 
further analysis.  But the ``refseq_protein`` database, which just contains a 
single reference sequence for each species (and no plasmids), doesn't have 
these problems in the first place.  In addition, I can have more confidence 
that the ``refseq_protein`` sequences are real proteins and don't have errors.
