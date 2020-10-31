# CONSULT
Accurate contamination removal using locality-sensitive hashing

CONSULT is the tool for contamination removal from genomic sequencing reads using exclusing read filtering approach. Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of the reference dataset. It supports the inclusion of approximately 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described design of the algorithm and software architecture will be available online shortly. <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

Summary data tables and scripts that we used during testing are available at [https://github.com/noraracht/lsh_scripts].
 <!--  - Raw data are deposited in -->


System Requirements
------------

**Disk space:** Construction of a CONSULT database requires approximately 120 GB of disk space. Exact footprint depends on a number of *k*-mers in a reference set. For instance, the size of the databases that we built during testing with default settings varied between 107G to 120GB. In GTDB database 62 GB of disk space was used to store encodings, 56 GB was taken by lookup table, 2 G was utilized by tag array and metadata. 

**Memory:** CONSULT is designed to operate on a machine with 128 GB of RAM. To run, it requires enough free memory to hold the database in RAM. we note that durign datatase construction the user will need slightly more than that in RAM to occomodate intermediary processes.

 
Installation
------------

CONSULT is a command-line tool implemented in c++11. It runs [Jellyfish][2] and [Mash][3] internally to efficiently compute k-mer profile of genome-skims and their intersection, and estimates the genomic distances by correcting for the effect of low coverage and sequencing error. Skmer also depends on [seqtk][5] for some FASTQ/A processings. 


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

