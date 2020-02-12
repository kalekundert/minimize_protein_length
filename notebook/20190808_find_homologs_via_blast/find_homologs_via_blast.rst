***********************
Find homologs via BLAST
***********************

The first step in this pipeline is to identify homologs of the protein of 
interest.  The most prominent ways to do this are BLAST and HMMer.  I decided 
to start with BLAST simply because it's more accessible from Biopython.

Considerations
==============

Remote vs. Local
----------------
BLAST searches can either be performed using the NCBI web server, or using a 
local copy of the sequence databases.  I decided to use the web server:

- This makes it easier to experiment with querying against different databases, 
  because I don't need to download each database I want to search against.  

- Additionally, if I put this project down for a while and come back to it 
  later, I won't need to re-download the databases to get the most up-to-date 
  sequences.

- Biopython has an API that makes BLAST queries using the web interface, so 
  it's easy to make these searches reproducible.  (Biopython has an API for 
  local BLAST searches too, but what I'm trying to avoid here is requiring the 
  user to visit the NCBI website and have to run queries through the web 
  interface.)

- Remote queries are significantly lower: 5-45 minutes as compared to seconds 
  (probably) for local queries.  But this is mitigated by the fact that (i) I 
  don't do many of these searches, (ii) I can do any number of searches 
  concurrently because none of my CPU/memory are used, (iii) the aforementioned 
  advantages of remote queries.

PSI-BLAST
---------
PSI-BLAST is an iterated version of BLAST that can be used to find more distant 
homologs.  My current impression is that this is not necessary for this 
application, as we are down-weighting distant homologs anyways.

Accession number vs. sequence
-----------------------------
Initially I used an accession number to specify the sequence to find homologs 
for.  But I switched to providing a protein sequence as input once I realized 
that I'd need to validate whether the sequence corresponding to an accession 
number actually matched the sequence being used in the assay.  I got the 
protein sequences by translating DNA sequences for each protein.  I got these 
DNA sequence from Xiao, and stored them as ``*.fasta`` files in the 
corresponding workspace directories.

Database
--------
There are a number of NCBI databases that BLAST can search.  I chose to use 
``refseq_protein`` instead of ``nr`` because ``nr`` has a lot of synthetic 
plasmids and nearly identical sequences from a handful of species (presumably 
the ones that have been sequenced multiple times).  These sequences don't add 
any information about which residues we can delete, so if I were to use ``nr``, 
I'd want to filter them out before doing any further analysis.  But the 
``refseq_protein`` database, which just contains a single reference sequence 
for each species (and no plasmids), doesn't have these problems in the first 
place.  In addition, I can have more confidence that the ``refseq_protein`` 
sequences are real proteins and don't have errors.


Methods
=======
::

   $ minp 02 sp_cas9
   Found 4891 homologs.

