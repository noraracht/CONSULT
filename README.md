# CONSULT
<!-- Accurate contamination removal using locality-sensitive hashing-->

CONSULT is the tool for contamination removal from genomic sequencing reads using exclusion read filtering approach. Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of the reference dataset. It supports the inclusion of approximately 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described design of the algorithm and software architecture will be available online shortly. <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

Summary data tables and scripts that we used during testing are available at https://github.com/noraracht/lsh_scripts.
 <!--  - Raw data are deposited in -->


System Requirements
------------

- **Disk space:** Construction of CONSULT database requires approximately 120 GB of disk space. Exact footprint depends on a number of *k*-mers in a reference set. For instance, the size of the databases that we built during testing with default settings varied between 107 GB to 120 GB. In GTDB database 62 GB of disk space was used to store encodings, 56 GB was taken by lookup table, 2 GB was utilized by tag array and metadata. 

- **Memory:** CONSULT is designed to operate on a machine with 128 GB of RAM. To run, it requires enough free memory to hold the entire database in RAM. We note that during datatase construction the user will need slightly more than that in RAM to accomodate intermediary processes.

 
Installation
------------
<!-- may be add download links at some point -->
- There is no installation required to run CONSULT.

- CONSULT is a command-line tool implemented in C++11 with some x86 assembly code. Database reading and query search are parallelized using [OpenMP](https://www.openmp.org). Core programs for map construction and query search need to be compiled using g++ 
```
    g++ main_map_v17.1.cpp -std=c++11 -O3 -o main_map
    g++  main_search_v17.4.cpp -std=c++11 -fopenmp -O3 -o main_search
```    

Using CONSULT
------------

<!--Change to the CONSULT working directory and run the scripts below. -->
 ### Database construction
To construct standard reference database, you can use the following command:
```
 ./main_map -i $INPUT_FASTA_FILE -o $DBNAME
```  
###### Input:
Input is supposed to be a FASTA file formatted as shown below. Specifically CONSULT is designed to accept [Jellyfish](http://www.genome.umd.edu/jellyfish.html) output files that represent a list of **32 bp** *k*-mers associated with their counts. See Data Preprocessing section for details on how to generate these files. Note, however, CONSULT does not use the count values and the only relevant information is the sequence itself. Jellyfish output is pseudo-randomly ordered, and thus further randomization is not needed. The sequences may be repeated. Duplicate entries will not be included in a database but since they are read and processed their presence will likely increase library construction time.

Example FASTA:
```
> FASTA sequence 1 label
AGACGAGCTTCTTCATCTAAAATGAATTCTCC
> FASTA sequence 2 label
CCAGCTGCATTGATGGTGGGAGTTGGTAAAGG
> FASTA sequence 3 label
GGACCTTGATTTTGACAAGATAGTCGGTAGAC
> FASTA sequence 4 label
ACCACATTTTATACATCGTAAGACAAGCGGCT
```

###### Output: 
Replace "$DBNAME" above with your preferred database name. Reference library will be created in the same directory where script is ran. If this working directory already contains a database with the same name software will throw an exception. This feature is included to prevent existing databases from being overwritten.

 ### Query search
To query a set of sequences against reference use the CONSULT command:
```
 ./main_search -i $DBNAME -c 0 -t 24 -q $QUERY_FOLDER
``` 
###### Input: 
The files containing query sequences to be classified should be located in $QUERY_FOLDER and be in a FASTQ format (one uncompressed .fq/.fastq file per each sample). FASTA format is not supported at the moment. Note, if you need to query FASTA files you can convert .fasta/.fa to .fastq/.fq using [seqtk](https://github.com/lh3/seqtk) ```"seqtk seq -F CHAR"``` command which attaches fake quality scores to the sequences. Quality factors are not being utilized by CONSULT but FASTQ labels will be used to identify the sequences in the output file.

Example FASTQ:
```
@ FASTQ sequence 1 label
CATCGAGCAGCTATGCAGCTACGAGT
+
-$#'%-#.&)%#)"".)--'*()$)%
@ FASTQ sequence 2 label
TACTGCTGATATTCAGCTCACACC
+
,*#%+#&*$-#,''+*)'&.,).,
```

###### Output: 
CONSULT is designed for filtering out contaminants from sequencing reads so its output is a FASTQ file that contains **unclassified** (clean) reads and their corresponding sequence IDs, obtained from the input FASTQ headers. Files are stored into working directory where software is ran. Every sample retains its original file name prefixed with *"ucseq_"*. 
<!--Log output is sent to standard output by default. -->

**CONSULT program arguments are:**

- -i - name of the reference database

- -c - the highest number of *k*-mers that is required to still keep sequencing read unclassified. For instance, if at least one *k*-mer match is enough to classify a read (default setting mentioned in a paper), "c" should be set to 0 in a software.  If at least two *k*-mer matches are required to call entire read a match, "c" should be set to 1. By default the value of "c" is 0.

- -t - number of threads

- -q - name of the folder where queries are located


Data Preprocessing
------------

We suggest the following workflow to obtain *k*-mer list file to construct CONSULT database from multiple assembly references:

**1. To combine assembly references** into single file:
```
cat /path/to/folder/*.fna > combined.fna
```


**2. Extraction of canonical 35 bp *k*-mers** from fasta genomic reference was performed with [Jellyfish](http://www.genome.umd.edu/jellyfish.html). 
 - To compute *k*-mer profile:
```
jellyfish count -m 35 -s 100M -t 24 -C combined.fna -o counts.jf
```
 - To output a list of *k*-mers associated with their counts:
 ```
jellyfish dump counts.jf > 35bp_kmer_lst.fa
```

**3. Minimization was performed using custom C++11 script.**  The script accepts as an input [Jellyfish](http://www.genome.umd.edu/jellyfish.html) fasta file containing 35 bp canonical *k*-mers extracted from reference and outputs their 32 bp minimizers in fasta format. 
- To compile:
```
g++ minimization_v3.0.cpp -std=c++11 -o main_minimization
```
- To run:
```./main_minimization -i 35bp_kmer_lst.fa -o 32bp_minzer_lst.fa```

**4. Extraction of canonical 32 bp *k*-mers**. 
- To compute *k*-mer profile and 
```
jellyfish count -m 32 -s 100M -t 24 -C 32bp_minzer_lst.fa -o counts.jf
```
 - To output a list of *k*-mers with their counts
 ```
jellyfish dump counts.jf > 32bp_kmer_lst.fa
```
 
We note, in our testing we used minimization technique to reduce *k*-mer count in original dataset. Alternatively if dataset is small and minimization is not needed the user can use last command directly to obtain a list of all 32 bp *k*-mers and utilize it as as an input into CONSULT software.



Quick start
------------
Share commands with toy examples

