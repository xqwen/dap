## TORUS - QTL Discovery utilizing Genomic Annotations 


``torus`` is a free software package that implements a computational procedure for discovering molecular QTLs incorporating genomic annotations. 


## Compilation

The following C/C++ libraries are required for compiling the source code

* GNU [GSL](http://www.gnu.org/software/gsl/) library
* [Zlib](http://zlib.net/) library
* [Boost](http://www.boost.org/) C++ library


## Command-line Syntax 

``torus -d data_file.gz -est|-qtl|-dump_prior [-annot annotation_file.gz] [-smap snp_map.gz] [-gmap gene_map.gz]``  

## Command line options

### 1. Required input files

* ``-d data_file.gz``: the data_file should be compressed in gzip format and  contain summary statistics from single SNP association analysis. Currently, two formats are supported: 1) gzipped MatrixEQTL output single SNP association result; 2) pre-computed Bayes factors based on single SNP association tests. 3) Summary-level Z-statistics from single SNP association tests. If 2) or 3) is used, it is important to specific ``--load_bf`` or ``--load_zval`` in the command line, respectively.



### 2. Optional input files

* ``-smap snp_map.gz``: gzipped SNP position files
* ``-gmap gene_map.gz``: gzipped Gene position (transcription start site) file
With both files specified, torus will control for SNP distance to transcript start site (TSS)


### 3. Task options


* ``-est``: output the estimates of enrichment parameters and their confidence intervals
* ``-qtl``: perform Bayesian FDR control, and output the result
* ``-dump_prior``: compute SNP-level priors using the estimated enrichment estimates for each locus

