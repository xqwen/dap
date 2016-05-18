# TORUS - QTL Discovery utilizing Genomic Annotations 


``torus`` is a free software package that implements a computational procedure for discovering molecular QTLs incorporating genomic annotations. 


## Compilation

The following C/C++ libraries are required for compiling the source code

* GNU GSL library
* Zlib library
* Boost C++ library


## Command line options

### 1. Required input files

* ``-d data_file``: the data_file should be compressed in gzip format and  contain summary statistics from single SNP association analysis. Currently, two formats are supported: 1) gzipped MatrixEQTL output single SNP association result;
2) pre-computed Bayes factors based on single SNP association tests. If 2) is used, it is important to specific ``--load_bf`` in the command line.



### 2. Optional input files

* ``-smap snp_map``: gzipped SNP position files
* ``-gmap gene_map``: gzipped Gene position (transcription start site) file
With both files specified, torus will control for SNP distance to transcript start site (TSS)


### 3. Task options


* ``-est``: output the estimates of enrichment parameters and their confidence intervals
* ``-qtl``: perform Bayesian FDR control, and output the result
* ``-dump_prior``: compute SNP-level priors using the estimated enrichment estimates for each locus


## Tutorials and Sample files

We provide two real data tutorials to illustrate the usage of ```torus``` for QTL discovery. They eQTL data from the GEUVADIS project and the GTEx liver tissue. The tutorials can be found in the [```examples```] (https://github.com/xqwen/torus/tree/master/examples/) directory. 


## Reference and citation

* Wen, X. Effective QTL Discovery Incorporating Genomic Annotations. *bioRxiv* doi:10.1101/032003.


