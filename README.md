# CONSULT
<!-- Accurate contamination removal using locality-sensitive hashing-->

CONSULT is a tool for contamination removal from genomic sequencing reads. Relying on locality-sensitive hashing, CONSULT extracts *k*-mers from a query set and tests whether they fall within a user-specified hamming distance of *k*-mers in the reference dataset. It supports the inclusion of approximately 8 billion *k*-mers in its reference library, accommodating datasets with tens of thousands of microbial species.

The paper where we have described design of the algorithm and software architecture will be available online shortly. <!-- (open access): -->
<!--  - [paper reference and doi][1] -->

* Summary data tables and scripts that we used during testing are available at https://github.com/noraracht/lsh_scripts.

* Raw data are deposited in https://github.com/noraracht/lsh_raw_data.


System Requirements
------------

- **Disk space:** Construction of CONSULT database requires approximately 120GB of disk space. Exact footprint depends on a number of *k*-mers in a reference set. For instance, the size of the 
