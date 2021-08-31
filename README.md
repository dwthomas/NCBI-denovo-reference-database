# MiFish-Reference-Database
Database of MiFish sequences from NCBI nucleotide database

The goal of this project is to mine MiFish seqeunces from the NCBI, and produce a database of sequences appropriate for use with Qiime2.
A key goal is to implement this such that the database can be incrementally updated as new sequences are added to the NCBI.

## Protocol
1. Query the ncbi nucleotide database, and download the list of matching accessions. I use the query "(((12S OR MiFish) OR mitochondrion) NOT chromosome) NOT shotgun".
2. For each accession fetch the sequence and taxonomic information, which can both be found by downloading the xml results from the nucleotide database with efetch.  Taxonomy is reformatted in the Qiime style, with a semicolon separated list of the lineage, followed by the scientific name, followed by a common name if available.  There are a handful of records without proper taxonomic information in the nucleotide database entries, those are queried from the taxonomy database.  The handful not resolved by that are excluded, they appear to all be short primer sequences and not useful to include.
3. The fetched sequences are length filtered to exclude any sequences too short to be useful, I filter at 100 bp.
4. We want to pull out the MiFish region of these references, there are plenty of mitochondrial genomes at this point and we want to subsequence those to only contain the MiFish amplified section.  To facilitate this I assembled the unique MiFish sequences found across our projects, these sequences represent a broad swath of fish life.  These sequences were blasted against each potential database sequence, the best hit selected and the aligned region extracted from the sequence.  Alignments longer than 140 bp were retained.
5. Sequences for each taxonomy string are dereplicated, this mainly massively reduces the number of identical sequences for a few highly studied organisms (ie Humans).  If the same sequence had multiple possible taxonomic strings those are retained.

## Implementation
1,2 are implemented by download_ncbi.py
3,4,5 are implemented in Taxa clarify.py

## Some notes
Downloading the whole database takes days.
The Vsearch classifier performs global alignments, so for it to work well we need to try and extract the actual MiFish region, and not just include a bunch of whole mitochondrial genomes.
This should be used for a fast initial taxonomy, I strongly recommend manually assigning taxonomy when using MiFish, as there are many examples of sequences that are ambiguous, but can be resolved with prior local knowledge of the species you would expect to see in your sample. 
