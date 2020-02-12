******************************
Find homologs in viral genomes
******************************

I had the thought that it might be beneficial to specifically look for viral 
homologs of the protein to be shortened.  The reason is that viral proteins 
(unlike most other proteins) are optimized to minimize length, which is what 
we're trying to do.  Of course, this might also dramatically limit the number 
of homologs found, limiting the utility of this approach.

Methods
=======

.. protocol::

   ::

      $ mkdir sp_cas9/blast_viral
      $ cd !$
      $ ln ../blast_refseq_5000/query.fasta .

   ``query.fasta`` is the protein sequence of the Cas9 we're using for the 
   assay.

   - Go to `NCBI Virus 
     <https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/sequence>`_.

   - Click on "Protein Search"

   - Copy the contents of ``query.fasta`` into the search form.  Remove the 
     ``*`` representing the stop codon.  Click "Start search".

   - When the search completes, download the results.

      - Click "Download"
      - Select data type: Protein
      - Select records: Download all records
      - Select FASTA definition line: Use default

   - Save the results as ``homologs.fasta`` in the ``blast_viral`` directory.

BLAST turned up 53 proteins with homology to Cas9 in viral genomes.  All of the 
results were poor.  Most aligned with about 5% of Cas9, and had only 20% 
identity in that region.  Furthermore, most of the results were on pol genes 
(e.g. reverse transcriptases), which I wouldn't expect to have any meaningful 
homology to Cas9.

I was more hoping to find phage genomes that had repurposed CRISPR proteins for 
something.  But it occurs to me that this database might only have human 
viruses (or at least, only eukaryotic viruses).

I'm not going to take this any further, because I don't think the results have 
any promise.

