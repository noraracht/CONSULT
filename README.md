# CONSULT
Accurate contamination removal using locality-sensitive hashing

CONSULT is the tool for contamination removal from genomic sequencing reads using exclusion read filtering approach. Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of the reference dataset. It supports the inclusion of approximately 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described design of the algorithm and software architecture will be available online shortly. <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

Summary data tables and scripts that we used during testing are available at https://github.com/noraracht/lsh_scripts.
 <!--  - Raw data are deposited in -->


System Requirements
------------

#### Disk space
Construction of CONSULT database requires approximately 120 GB of disk space. Exact footprint depends on a number of *k*-mers in a reference set. For instance, the size of the databases that we built during testing with default settings varied between 107 GB to 120 GB. In GTDB database 62 GB of disk space was used to store encodings, 56 GB was taken by lookup table, 2 GB was utilized by tag array and metadata. 

 #### Memory
CONSULT is designed to operate on a machine with 128 GB of RAM. To run, it requires enough free memory to hold the entire database in RAM. We note that during datatase construction the user will need slightly more than that in RAM to accomodate intermediary processes.

 
Installation
------------

CONSULT is a command-line tool implemented in C++11 with some x86 assembly code. Database reading and query search are parallelized. Multithreading is handled using [OpenMP](https://www.openmp.org). 
  
Core programs for map construction and query search need to be compiled using g++ 
```
    g++ main_map_v17.1.cpp -std=c++11 -O3 -o main_map
    g++  main_search_v17.4.cpp -std=c++11 -fopenmp -O3 -o main_search
```    

Using CONSULT
------------

<!--Change to the CONSULT working directory and run the scripts below. -->
#### Database construction
To construct standard reference database, you can use the following command:
```
 ./main_map -i $INPUT_FASTA_FILE -o $DBNAME
```  
Replace "$DBNAME" above with your preferred database name. Reference library will be created in the same directory where script is ran. If this working directory already contains a database with the same name software will throw an exception. This feature is included to prevent existing database from being overwritten.

#### Query search
To query a set of sequences against reference use the CONSULT command:
```
 ./main_search -i $DBNAME -c 0 -t 24 -q $QUERY_FOLDER
``` 
Output will be sent to standard output by default **!!!!- not correct!!**. The files containing the sequences to be classified should be located in $QUERY_FOLDER and be in a FASTQ format (one uncompressed .fq/.fastq file per each sample). FASTA format is not supported at the moment. However, if you need to query FASTA files you can convert .fasta/.fa to .fastq/.fq using ** [samtools]() !!! check** which attached dummy quality score to the sequences.

**CONSULT program arguments are:**

- -i - name of the reference database

- -c - the highest number of *k*-mers that is required to still keep sequencing read unclassified. For instance, if at least one *k*-mer match is enough to classify a read, "c" should be set to 0.  If at least two *k*-mer matches are required to call entire read a match, "c" should be set to 1. Default setting for c in 0.

- -t - number of threads

- -q - name of the folder where queries are located


Data Preprocessing
------------
CONSULT accepts as an input *k*-mer output file from [Jellyfish](http://www.genome.umd.edu/jellyfish.html)
1. Cat fna files
2. jellyfish
3. minimization
4. jellyfish again (to remove canonical \kmers so they don't talke up db space)

<!--It runs [Jellyfish][2] and [Mash][3] internally to efficiently compute k-mer profile of genome-skims and their intersection, and estimates the genomic distances by correcting for the effect of low coverage and sequencing error. Skmer also depends on [seqtk][5] for some FASTQ/A processings. -->

```
g++-9 minimization_v3.0.cpp -std=c++11 -o main_minimization
```


Quick start
------------

