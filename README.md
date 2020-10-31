# CONSULT
Accurate contamination removal using locality-sensitive hashing

CONSULT is the tool for contamination removal from genomic requencing reads using exclusing read filtering approach. The paper where we have described design of the algorithm and software architecture will be available online shortly. <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

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

