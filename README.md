# CONSULT
<!-- Accurate contamination removal using locality-sensitive hashing-->

CONSULT is a tool for contamination removal from genomic sequencing reads. Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of *k*-mers in the reference dataset. It supports the inclusion of approximately 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described design of the algorithm and software architecture are now available on online (https://www.biorxiv.org/content/10.1101/2021.03.18.436035v1). <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

* Summary data tables and scripts that we used during testing are available at https://github.com/noraracht/lsh_scripts.

* Raw data are deposited in https://github.com/noraracht/lsh_raw_data.

* Our custom CONSULT libraries constructed using different genomic reference sets:
    - [GTDB](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_gtdb.tar.gz)
    - [TOL](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_ToL.tar.gz)
    - [Bacterial/Archaeal Kraken](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_kraken.tar.gz)
    - [Mitochondrial CONSULT database](https://tera-trees.com/data/consult/v1.0.0/consult_mito_k32_p3l2clmn7_K15_tag2_v171.tar.gz)
  
  <!--  - [GTDB](https://drive.google.com/file/d/1MQJAXmZiTurumlZpvNoMLB0tKWGM_VE4/view?usp=sharing)-->
  <!--    - [TOL](https://drive.google.com/file/d/1sA9HFjWoU2jZ2vjd98pHVDEFRzOKMImk/view?usp=sharing)-->
  <!--    - [Bacterial/Archaeal Kraken](https://drive.google.com/file/d/1jeZB6b6aXl06BpPPsjM8oQA4xingJ1Dq/view?usp=sharing)-->
  <!--   - [Mitochondrial CONSULT database](https://drive.google.com/file/d/1mFD3dYFrJKqUkWlkRHbrQt-6eG-_K5vI/view?usp=sharing)-->

At the moment when using our libraries, library name needs to stay unchanged since library files are prefixed with library name. This will change in future releases to allow for more flexibility.


System Requirements
------------

- **Disk space:** Construction of CONSULT database requires approximately 120GB of disk space. Exact footprint depends on a number of *k*-mers in a reference set. The size of the three main databases that we built for testing with default settings varied between 107GB to 120GB. For example, for the GTDB database, 62GB of disk space was used to store encodings, 56GB was taken by lookup table, and 2GB was utilized by tag array and metadata. 

- **Memory:** CONSULT is designed to operate on a machine with 128GB of RAM. To run, it requires enough free memory to hold the entire database in RAM. We note that during datatase construction the user will need slightly more RAM than 128GB to accomodate intermediary processes. Once the database is built, 128GB should be sufficient. 

- **Dependencies:** CONSULT is a command-line tool implemented in C++11 with some x86 assembly code. Database reading and query search are parallelized using [OpenMP](https://www.openmp.org). 
    - Core programs for map construction and query search need to be compiled using somewhat recent version of g++ that will support C++11. For our tests we compiled with versions 4.8.5 and 7.2.0. 
    - The database construction uses some external tools such as [Jellyfish](http://www.genome.umd.edu/jellyfish.html).

 
Installation
------------
<!-- may be add download links at some point -->

1. Download using one of two approaches:
    - You can obtain the [zip file](https://github.com/noraracht/CONSULT/archive/main.zip) and extract the contents to a folder of your choice. Then, proceed to compilation.
    - Alternatively, you can clone the [github repository](https://github.com/noraracht/CONSULT.git) and continue with compilation.

<!--- CONSULT is a command-line tool implemented in C++11 with some x86 assembly code. Database reading and query search are parallelized using [OpenMP](https://www.openmp.org). Core programs for map construction and query search need to be compiled using g++ -->
2. To compile, go to the directory where core programs for map construction and query search are located and run:
```
    g++ main_map.cpp -std=c++11 -O3 -o main_map
    g++ main_search.cpp -std=c++11 -fopenmp -O3 -o main_search
```    

Execution
------------

<!--Change to the CONSULT working directory and run the scripts below. -->
 ### Database construction
To construct standard reference database, go to the place where scripts were compiled <!--software was builf--> and use the following command:
```
 ./main_map -i $INPUT_FASTA_FILE -o $DBNAME
```  
###### Input:
Input is supposed to be a FASTA file formatted as shown below. Specifically CONSULT is designed to accept [Jellyfish](http://www.genome.umd.edu/jellyfish.html) output files that represent a list of **32 bp** *k*-mers associated with their counts. We tested with [Jellyfish](http://www.genome.umd.edu/jellyfish.html) 2.3.0. See Data preprocessing section for details on how to generate the input file.  Note that CONSULT does not use the count values and the only relevant information is the sequence itself. Jellyfish output is pseudo-randomly ordered, and thus further randomization is not needed. The sequences should not be repeated.

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
Replace "$DBNAME" above with your preferred database name. Reference library will be created in the same directory where the script is run. If this working directory already contains a database with the same name, the software will throw an exception. This feature is included to prevent existing databases from being overwritten.

 ### Query search
To query a set of sequences against reference, go to the directory where binaries are and execute the CONSULT command:
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
CONSULT is designed for filtering out contaminants from sequencing reads so its output is a FASTQ file that contains **unclassified** (clean) reads and their corresponding sequence IDs, obtained from the input FASTQ headers. Files are stored into working directory where software is run. Every sample retains its original file name prefixed with *"ucseq_"*. 
<!--Log output is sent to standard output by default. -->

**CONSULT program arguments are:**

- -i - name of the reference database

- -c - the lowest number of *k*-mers that is required to call sequencing read classified. For instance, if at least one *k*-mer match is enough to classify a read (default setting mentioned in a paper), "c" should be set to 1 in a software.  If at least two *k*-mer matches are required to call entire read a match, "c" should be set to 2. By default the value of "c" is 1.

- -t - number of threads

- -q - name of the folder where queries are located


Data preprocessing
------------

We suggest the following workflow to obtain *k*-mer list file to construct CONSULT database from multiple assembly references.

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
g++ minimization.cpp -std=c++11 -o main_minimization
```
- To run:
```
./main_minimization -i 35bp_kmer_lst.fa -o 32bp_minzer_lst.fa
```

**4. Extraction of canonical 32 bp *k*-mers** allowed to further reduce *k*-mer count and remove duplicate sequences. 
- To compute *k*-mer profile:
```
jellyfish count -m 32 -s 100M -t 24 -C 32bp_minzer_lst.fa -o counts.jf
```
- To output a list of *k*-mers with their counts:
``` 
jellyfish dump counts.jf > 32bp_kmer_lst.fa
```
 
We note, in our testing we used minimization technique to reduce *k*-mer count in original dataset. Alternatively, if dataset is small and minimization is not needed the user can use last command directly to obtain a list of all 32 bp *k*-mers and utilize it as an input into CONSULT software.


Quick start
------------
To allow users to get familiar with the software and test their installations we have provided a few example files. After downloading and compiling the software, go to the directory where CONSULT was installed and run:

### Minimization testing
```
./main_minimization -i k35C_bef_mininimization.fa -o k32C_af_mininimization.fa
```
Sample [k35C_bef_mininimization.fa](https://github.com/noraracht/CONSULT/blob/main/k35C_bef_mininimization.fa) contains 600000 35 bp *k*-mer sequences and takes approximately 10s to run locally. Minimized sequences are stored in [k32C_af_mininimization.fa](https://github.com/noraracht/CONSULT/blob/main/k32C_af_mininimization.fa).

### Library construction testing
```
./main_map -i k32C_af_mininimization.fa -o G000307305_nbr_map
```
This step takes about 7 min to complete. The constructed database uses ~60GB of disk space.

### Query search testing
```
./main_search -i G000307305_nbr_map -c 1 -t 1 -q query_set
```
Sample query [G000307305.fq](https://github.com/noraracht/CONSULT/blob/main/query_set/G000307305.fq) contains 66667 genomic reads. Classification running time is ~2 min with 1 thread. Approximately 38000 reads (in our case 38440) from the query should match to the database. Since library construction involves randomization exact number of matched sequences will be slightly different. Unclassified reads are stored in [ucseq_G000307305.fq](https://github.com/noraracht/CONSULT/blob/main/ucseq_G000307305.fq).


On a side note, during search for every query sample CONSULT outputs to the standard output: sample name, sample line count and number of matched reads. Since each entry in a FASTQ file consists of 4 lines, total number of entry sequences as well as percentage of matched reads can be easily computed using these values.


 <!--We have tested CONSULT only on Linux and MAC. To report issues please use-->
