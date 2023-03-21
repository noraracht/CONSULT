# CONSULT

CONSULT is a tool for contamination removal from genomic sequencing reads.
Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of *k*-mers in the reference dataset.
It supports the inclusion of approximately 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described the design of the algorithm and software architecture is now available in the [publication](https://academic.oup.com/nargab/article/3/3/lqab071/6342218)).

* Summary data tables and scripts that we used during testing are available at [GitHub](https://github.com/noraracht/lsh_scripts).

* Raw data are deposited in [this repository](https://github.com/noraracht/lsh_raw_data).

* Our custom CONSULT libraries are constructed using different genomic reference sets:
    - [GTDB](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_gtdb.tar.gz)
    - [TOL](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_ToL.tar.gz)
    - [Bacterial/Archaeal Kraken](https://tera-trees.com/data/consult/v1.0.0/all_nbrhood_kmers_k32_p3l2clmn7_K15-map2-171_kraken.tar.gz)
    - [Mitochondrial CONSULT database](https://tera-trees.com/data/consult/v1.0.0/consult_mito_k32_p3l2clmn7_K15_tag2_v171.tar.gz)
    - [Plastid CONSULT database](https://tera-trees.com/data/consult/v1.0.0/consult_plastid_k32_p3l2clmn7_K15_tag2_v171.tar.gz)

 *  Alternatively our custom CONSULT libraries can be downloaded from the cloud drive:
    - [GTDB](https://drive.google.com/file/d/1MQJAXmZiTurumlZpvNoMLB0tKWGM_VE4/view?usp=sharing)
    - [TOL](https://drive.google.com/file/d/1sA9HFjWoU2jZ2vjd98pHVDEFRzOKMImk/view?usp=sharing)
    - [Bacterial/Archaeal Kraken](https://drive.google.com/file/d/1jeZB6b6aXl06BpPPsjM8oQA4xingJ1Dq/view?usp=sharing)
    - [Mitochondrial CONSULT database](https://drive.google.com/file/d/1mFD3dYFrJKqUkWlkRHbrQt-6eG-_K5vI/view?usp=sharing)


System requirements
------------

- **Disk space:** Exact footprint depends on a number of *k*-mers in a reference set, and the parameters used.
The size of the three main databases that we built for testing without using the default heuristic, using the parameter configuration described in the [publication](https://academic.oup.com/nargab/article/3/3/lqab071/6342218), varied between 107GB to 120GB.
For example, for the GTDB database, 62GB of disk space was used to store encodings, 56GB was taken by lookup table, and 2GB was utilized by tag array and metadata.
CONSULT will report the estimated library size during the mapping.
For instance, $2^{30}$ *k*-mers will result in the constructed library being approximately 17GB, if the default heuristic is used.

- **Memory:** CONSULT requires enough free memory to hold the entire database in RAM.
Using the default heuristic, some examples of approximately required memory with respect to the number of *k*-mers in a reference set can be listed as below;
    * $2^25$ *k*-mers $\rightarrow$ 600MB,
    * $2^{26}$ *k*-mers $\rightarrow$ 1.2GB,
    * $2^{28}$ *k*-mers $\rightarrow$ 4.3GB,
    * $2^{30}$ *k*-mers $\rightarrow$ 17.1GB,
    * $2^{32}$ *k*-mers $\rightarrow$ 68.7GB.

We note that during library construction the user will need slightly more RAM than the given values to accommodate intermediary processes.
Once the database is built, these values should be sufficient.

- **Dependencies:** CONSULT is a command-line tool implemented in C++11 with some x86 assembly code. Database reading and query search are parallelized using [OpenMP](https://www.openmp.org).
    - We complied map construction and query search with a g++ that supports C++11 (required).
    For our tests, we have compiled versions 4.8.5 and 7.2.0, both of which work.
    - The **database construction** uses some external tools such as [Jellyfish](http://www.genome.umd.edu/jellyfish.html).


Installation
------------

1. Download using one of two approaches:
    - You can obtain the [zip file](https://github.com/noraracht/CONSULT/archive/main.zip) and extract the contents to a folder of your choice.
    Then, proceed to compilation.
    - Alternatively, you can clone the [github repository](https://github.com/noraracht/CONSULT.git) and continue with compilation.

2. To compile, go to the directory where core programs for map construction and query search are located and run the below commands.
    * You can use ``make`` to compile CONSULT.
    ```bash
    make all # for all components of CONSULT
    # OR
    make minimize # for the minimization script
    make map # for consult_map to construct a library
    make search # for consult_search to make queries
    ```
    * Alternatively, you can run ``g++`` directly.
    ```bash
    g++ minimize.cpp -std=c++11 -o minimize # for the minimization script
    g++ consult_map.cpp -std=c++11 -O3 -o consult_map  # for consult_map to construct a library
    g++ consult_search.cpp -std=c++11 -fopenmp -O3 -o consult_search # for consult_search to make queries
    ```

Constructing a library
-----

There are two steps:

1. Pre-processing your input genomes.
2. Construction of the CONSULT library.

### 1. Preprocessing
We suggest the following workflow to obtain the *k*-mer list file to construct the CONSULT database from multiple assembly references.

1. **To combine assembly references** into single file use ``cat`` as follows.
```bash
cat /path/to/folder/*.fna > combined.fna
```

2. **Extraction of canonical 35 bp *k*-mers** from FASTA genomic reference was performed with [Jellyfish](http://www.genome.umd.edu/jellyfish.html).
 - To compute *k*-mer profile run ``jellyfish count``.
```bash
jellyfish count -m 35 -s 100M -t 24 -C combined.fna -o counts.jf
```
 - To output a list of *k*-mers associated with their counts use ``jellyfish dump``.
 ```bash
jellyfish dump counts.jf > 35bp_kmer_lst.fa
```

3. **Minimization was performed using custom C++11 script.**
The script accepts as an input [Jellyfish](http://www.genome.umd.edu/jellyfish.html) fasta file containing 35 bp canonical *k*-mers extracted from reference and outputs their 32 bp minimizers in fasta format.
Run ``minimize`` as below.
```bash
./minimize -i 35bp_kmer_lst.fa -o 32bp_minzer_lst.fa
```

4. **Extraction of canonical 32 bp *k*-mers** allowed to further reduce *k*-mer count and remove duplicate sequences.
 - To compute *k*-mer profile run ``jellyfish count`` again.
```bash
jellyfish count -m 32 -s 100M -t 24 -C 32bp_minzer_lst.fa -o counts.jf
```
 - To export the list of *k*-mers and their counts use ``jellyfish dump``.
```bash
jellyfish dump counts.jf > 32bp_kmer_lst.fa
```

We note, in our testing we used the minimization technique to reduce the *k*-mer count in the original datasets.
Alternatively, if a dataset is small and minimization is not needed the user can use the last command directly to obtain a list of all 32 bp *k*-mers and utilize it as input into CONSULT software.

### 2. Mapping
To construct a standard reference database, go to the place where scripts were compiled and use the following command to run ``consult_map``.
```bash
 ./consult_map -i $INPUT_FASTA_FILE -o $DBNAME
```
###### Input:
Input is supposed to be a FASTA file formatted as shown below.
Specifically CONSULT is designed to accept [Jellyfish](http://www.genome.umd.edu/jellyfish.html) output files that represent a list of 32 bp *k*-mers associated with their counts.
We tested with [Jellyfish](http://www.genome.umd.edu/jellyfish.html) 2.3.0.
See the data preprocessing section for details on how to generate the input file.
Note that CONSULT does not use the count values and the only relevant information is the sequence itself.
Jellyfish output is pseudo-randomly ordered, and thus, further randomization is not needed.
The sequences should not be repeated.

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
Replace ``$DBNAME`` above with your preferred database name.
The reference library will be created in the same directory where the ``consult_map`` is run.
If this working directory already contains a database with the same name, the software will throw an exception.
This feature is included to prevent existing databases from being overwritten.

### A Quick Test

To allow users to get familiar with the software and test their installations we have provided a few example files. After downloading and compiling the software, go to the directory where CONSULT was installed and run:

### Minimization testing
```bash
cd example
./minimize -i k35C_bef_mininimization.fa -o k32C_af_mininimization.fa
# OR
make minimize # simply run this command in example/ directory
```

* Sample [k35C_bef_mininimization.fa](https://github.com/noraracht/CONSULT/blob/main/k35C_bef_mininimization.fa) contains 600000 35 bp *k*-mer sequences and takes ~3 seconds to run locally.
* Minimized sequences are stored in [k32C_af_mininimization.fa](https://github.com/noraracht/CONSULT/blob/main/k32C_af_mininimization.fa).

### Library construction testing
```bash
# in the example/ directory
./consult_map -p 3 -t 2 -i k32C_af_mininimization.fa -o G000307305_nbr_map
# OR
make map # simply test with make command in example/ directory
```

* This step takes about ~5 seconds to complete.
* The constructed database should use less than 500MB of disk space.
* ``-t`` is tag size in bits, and determines the number of partitions ( $2^t$ ).
* ``-p`` is the Hamming distance threshold for a match.


Searching against a library
------------
To query a set of sequences against a reference, go to the directory where binaries are and execute the CONSULT command:
```bash
 ./consult_search -i $DBNAME -q $QUERY_PATH -o $OUTPUT_DIRECTORY
```

###### Input:
The files containing query sequences to be classified should be located in ``$QUERY_PATH`` and be in a FASTQ format (one uncompressed ``.fq``/``.fastq`` file per each sample).
The path ``$QUERY_PATH`` can be a directory or a ``.fastq``.
If it is a directory, each query file in the directory will be queried against the library, and separate outputs will be generated for each.
FASTA format is not supported at the moment.
Note, if you need to query FASTA files you can convert ``.fasta``/``.fa`` to ``.fastq``/``.fq`` using [seqtk](https://github.com/lh3/seqtk) ``seqtk seq -F CHAR`` command which attaches fake quality scores to the sequences.
Quality factors are not being utilized by CONSULT but FASTQ labels will be used to identify the sequences in the output file.

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
CONSULT is designed for filtering out contaminants from sequencing reads.
So, its default output is a FASTQ file that contains **unclassified** reads and their corresponding sequence IDs, obtained from the input FASTQ headers.
Files are stored in the directory given in ``$OUTPUT_DIRECTORY``, and the default is where software is run.
Every sample retains its original file name prefixed with *"unclassified-seq_"*.

CONSULT also is able to generate a file that contains the **classified reads**, in the same format with **unclassified** described above, and the output file name will be prefixed with *"classified-seq_"*.
To make CONSULT behave this way, give the ``--classified-out`` flag.

Another output that CONSULT can generate is a tab-separated file, in which, each row is a read and the column values are the total number of *k*-mers in that read which matches with some *k*-mer in the reference library with Hamming distance $d$.
Each column corresponds to a distance value $d$.
CONSULT outputs this file, named file name prefixed with *"kmer-distances_*" if the ``--save-distances.`` flag is given.
The maximum distance value included in this file is determined by the $p$ value (Hamming distance threshold for a match) as a default ( $\lceil 1.5p \rceil$ ), but can be set to some other value with argument ``--maximum-distance``.

### A Quick Test of Query Search

* You have previously constructed a toy-example reference library called ``G000307305_nbr_map/``.
* Now we will test against that library using the query file called ``G000307305.fq``.

Use the following command to search ``G000307305.fq`` against ``G000307305_nbr_map``.

```bash
cd example
./consult_search  -c 1 --thread-count 1 -i G000307305_nbr_map -q G000307305.fq
# OR
make search # to use make for testing consult_search
```

* The simple query [G000307305.fq](https://github.com/noraracht/CONSULT/blob/main/query_set/G000307305.fq) contains 66667 genomic reads.
* Classification running time is ~5 seconds with 1 thread.
* Change ``--thread-count`` value to determine the number of threads, and benefit from parallelism.
Each read is searched in parallel, and the speed-up is almost linear with respect to the number of threads.
* ``-c`` is the minimum number of matching *k*-mers that is required to call a sequencing read classified.
* Approximately 38000 reads (in our case 38440) from the query should match the database.
* Unclassified reads are stored in ``unclassified-seq_G000307305.fq``.
* If flag ``--classified-out`` was given, classified reads would be stored in ``classified-seq_G000307305.fq``.
* If flag ``--save-distances`` was given, distance values would be stored in ``kmer-distances_G000307305.fq``.


Description of CONSULT arguments
--------------------------------

### ``minimize``

- ``-i`` or ``--input-fasta-file``: input ``.fasta`` file containing canonical *k*-mers, default length is 35.
- ``-o`` or ``--output-fasta-file``: ``.fasta`` file to output minimizers of the given canonical *k*-mers, default length is 32.

### ``consult_map``

- ``-i`` or ``--input-fasta-file``: input ``.fasta`` file to construct library.
- ``-o`` or ``--output-library-dir``: output path to the directory that will constitute the CONSULT library.
- ``-h`` or ``--number-of-positions``: number of randomly positioned bits to compute LSH.
- ``-t`` or ``--tag-size``: number of bits to be used as tag.
- ``-l`` or ``--number-of-tables``: number of tables, i.e., number of hash functions.
- ``-b`` or ``--column-per-tag``: number of columns per each tag partition, i.e., number of *k*-mers each encoding can map to.

### ``consult_search``
- ``-i`` or ``input-library-dir``: directory of the CONSULT library that will be used as the reference database.
- ``-o`` or ``output-result-dir``: directory in which classified reads, unclassified reads and matched *k*-mer counts will be saved.
- ``-q`` or ``--query-path``: the path to the query file, or to the directory containing query files.
- ``-c`` or ``--number-of-matches``: the minimum number of matched *k*-mers that is required to call sequencing read classified.
For instance, if at least one *k*-mer match is enough to classify a read (default setting mentioned in a paper), ``-c`` should be set to 1 in the software.
If at least two *k*-mer matches are required to call the entire read a match, ``-c`` should be set to 2.
Default value is 1.
- ``--thread-count``: number of threads to be used.
- ``--unclassified-out``: to output reads that are unclassified in a file with a name query file name prefixed with *"unclassified-seq_"*.
This is given by default.
- ``--classified-out``: to output reads that are classified in a file with a name query file name prefixed with *"classified-seq_"*.
- ``--save-distances``: to save number of matched *k*-mers in a tab-separated file where columns are the distances of corresponding counts.
File name is query file name prefixed with *"kmer-distances_"*.
- ``--maximum-distance``: maximum distance to be included as a column in the file containing *k*-mer match counts with respect varying Hamming distance values.
Note that, when the ``--maximum-distance`` value is too large compared to $p$, CONSULT does not necessarily aim to find such *k*-mers to compute column values.
This is because the library size and the corresponding theoretical guarantees only consider $p$ value.
So, if one would like to use a large ``--maximum-distance`` value, $p$ value should be increased proportionally.


Useful tips
------------

- **Large size queries:** When a larger plant or mammalian samples (>20 million reads) need to be queried there might be a need to speed up computation.
In such cases, we suggest splitting the query into smaller samples to process them independently on multiple nodes and combine outputs downstream.

- During a search for every query sample CONSULT outputs to the standard output: sample name, sample line count, and the number of matched reads.
Since each entry in a FASTQ file consists of 4 lines, the total number of entry sequences as well as the percentage of matched reads can be easily computed using these values.
