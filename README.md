# CONSULT
<!--  - Accurate contamination removal using locality-sensitive hashing -->

CONSULT is the tool for contamination removal from genomic sequencing reads using exclusing read filtering approach. Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of the reference dataset. It allows to include close to 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described design of the algorithm and software architecture will be available online shortly. <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

CONSULT extracts k-mers from sequencing reads and searches for an inexact match between k-mers of the query and a database. If match is found, That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.

  

CONSULT examines a constant number of k-mers in a a reference dataset and test whether they fall within some user-specified hamming distance fromthe reference se


CONSULT is a command-line tool implemented in c++11. It runs [Jellyfish][2] and [Mash][3] internally to efficiently compute k-mer profile of genome-skims and their intersection, and estimates the genomic distances by correcting for the effect of low coverage and sequencing error. Skmer also depends on [seqtk][5] for some FASTQ/A processings. 

System Requirements
------------

Blah


Installation
------------
**On 64-bit Linux and Mac OSX**, you can install  from bioconda channel using conda package manager. 
1. Install [Miniconda][4] (you can skip this if you already have either of Miniconda or Anaconda installed). 
2. Add the bioconda channel by running the following commands in your terminal (order matters):
```
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
```    


CONSULT Databases
------------

