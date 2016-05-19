## TORUS - QTL Discovery utilizing Genomic Annotations 


``torus`` implements the EM-DAP1 algorithm for enrichment analysis and QTL discovery.


## Compilation

The following C/C++ libraries are required for compiling the source code

* GNU [GSL](http://www.gnu.org/software/gsl/) library
* [Zlib](http://zlib.net/) library
* [Boost](http://www.boost.org/) C++ library


## Command-line Syntax 

``torus -d data_file.gz -est|-qtl|-dump_prior [-annot annotation_file.gz] [-smap snp_map.gz] [-gmap gene_map.gz]``  

## Command line options

### 1. Required input files

* ``-d data_file.gz``: the data_file should be compressed in gzip format and  contain summary statistics from single SNP association analysis. 

### 2. Task options

* ``-est``: perform joint enrichment analysis of annotations, output the estimates of enrichment parameters and their confidence intervals
* ``-qtl``: perform Bayesian FDR control for QTL discovery
* ``-dump_prior output_dir``: perform joint enrichment analysis of annotations and compute SNP-level priors using the estimated enrichment estimates for each SNP in each locus. The priors for each locus are saved in a single file in the ``output_dir`` and the file is ready for use by DAP to perform integrative fine-mapping analysis


### 3. Optional input files

* ``-annot annotation_file.gz``: the annotation_file should be compressed in gzip format and contain SNP annotation information 
* ``-smap snp_map.gz``: this option is designed for controlling SNP distance to transcript start site (TSS) in cis-eQTL mapping
* ``-gmap gene_map.gz``: this option is designed for controlling SNP distance to TSS in cis-eQTL mapping. The gene_map file contains TSS information of each target gene. With both ``-smap`` and ``-gmap`` specified, torus will automatically control SNP distance to TSS in enrichemnt analysis and QTL discovery. 

If none of the above three options are specified, torus simply treats each candidate SNP exchangeable and estimates only the baseline enrichemnt parameter (alpha_0).

## File formats


### 1. Input data format 

Currently, the following three input formats are supported:

* MatrixEQTL output of single-SNP association results, e.g., 
```
SNP	locus  beta	     t-stat	 p-value
rs1234	geneA  0.13	     1.28	 0.198    
...
``` 
A header is allowed (but not required) in the input file. This is the format directly output from the software [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/). This is the default input format recognized by torus.


* pre-computed Bayes factors based on single SNP association tests, e.g.,
```
SNP	locus	log10_BF
rs1234	geneA	0.234 
...
``` 
A header is allowed (but not required) in the input file. It is important to specify the command-line option ``-load_bf`` when this particular format is used in the input data file.
Using this  format, TORUS can deal with more complicated data structures, e.g., meta-analysis, multi-tissue eqtls, by pre-computing SNP-level Bayes factors from the softwares like [MeSH](https://github.com/xqwen/mesh) and [eQTLBMA](https://github.com/timflutre/eqtlbma). 



* Summary-level Z-statistics from single SNP association tests, e.g.,
```
SNP	locus	z-val
rs1234	geneA	1.28
...
```
A header is allowed (but not required) in the input file. It is	important to specify the command-line option ``-load_zval`` when this particular format is used in the input data file. This format is probably most convenient for analyzing GWAS data, for which single-SNP association z-values are typically available.


### 2. Annotation file format

The annotation requires a header to specify the name and the nature (categorical or continuous) of each SNP-level anntation. For example,
```
SNP   binding_d
chr1.51479  0
chr1.52058  2
chr1.52238  1
```
The first column  ("SNP") alwasys represents the SNP name. The following columns represent specific annotations. For categorical/discrete annotations, the annotation name should alwasy have a suffix "_d"; whereas for continuous annotations, the name should end with "_c".  In the above example, the annotation "binding" is a categorical variable with 3 categories (0, 1 and 2).  There is no restriction for the number of annotations for enrichment analysis.

 



