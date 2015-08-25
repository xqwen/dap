# DAP-1 embedded EM algorithm for enrichment analysis

This direcotory contains the source code for integrative enrichment analysis of QTL data using the DAP-1 algorithm. Although the exact functionality can be achieved using the [script](https://github.com/xqwen/dap/blob/master/experiments/enrichment/em_fmap.pl) with the ``dap`` executable, we provide this C++ implementation for much improved computational efficiency. The project is still under active development. 

## Compilation 

The following C/C++ libraries are required for compiling the source code

* GVU GSL library
* Zlib library
* Boost C++ library


## Command line options

### 1. Required input files

* ``-d data_file``: the data_file should be compressed in gzip format and  contain summary statistics from single SNP association analysis. Currently, two formats are supported: 1) gzipped MatrixeQTL output single SNP association result; 2) pre-computed Bayes factors based on single SNP association tests. If 2) is used, it is important to specific ``--load_bf`` in the command line.



### 2. Optional input files

* ``-smap snp_map``: gzipped SNP position files
* ``-gmap gene_map``: gzipped Gene position (transcription start site) file



## Sample files

The example to run ``em_dap1`` can be found in [here](https://github.com/xqwen/dap/tree/master/experiments/GEUVADIS)


    