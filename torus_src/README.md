## TORUS - QTL Discovery utilizing Genomic Annotations 


``torus`` implements the EM-DAP1 algorithm for enrichment analysis and QTL discovery.


## Compilation

The following C/C++ libraries are required for compiling the source code

* GNU [GSL](http://www.gnu.org/software/gsl/) library
* [Zlib](http://zlib.net/) library
* [Boost](http://www.boost.org/) C++ library


## Command-line Syntax 

``torus -d data_file.gz -est|-qtl|-dump_prior [-annot annotation_file.gz] [-smap snp_map.gz] [-gmap gene_map.gz] [-alpha fdr_level]``  

## Command-line Options

### 1. Required Input Files

* ``-d data_file.gz``: the data_file should be compressed in gzip format and  contain summary statistics from single SNP association analysis. 

### 2. Task Options

* ``-est``: perform joint enrichment analysis of annotations, output the estimates of enrichment parameters and their confidence intervals
* ``-qtl``: perform Bayesian FDR control for QTL discovery
* ``-dump_prior output_dir``: perform joint enrichment analysis of annotations and compute SNP-level priors using the estimated enrichment estimates for each SNP in each locus. The priors for each locus are saved in a single file in the ``output_dir`` and the file is ready for use by DAP to perform integrative fine-mapping analysis


### 3. Optional Input Files

* ``-annot annotation_file.gz``: the annotation_file should be compressed in gzip format and contain SNP annotation information 
* ``-smap snp_map.gz``: this option is designed for controlling SNP distance to transcript start site (TSS) in cis-eQTL mapping
* ``-gmap gene_map.gz``: this option is designed for controlling SNP distance to TSS in cis-eQTL mapping. The gene_map file contains TSS information of each target gene. With both ``-smap`` and ``-gmap`` specified, torus will automatically control SNP distance to TSS in enrichemnt analysis and QTL discovery. 

If none of the above three options are specified, torus simply treats each candidate SNP exchangeable and estimates only the baseline enrichemnt parameter (alpha_0).


### 4. Other options

* ``-alpha fdr_control_level``: pre-defined FDR control level. By default, it is set to 0.05

## Input File format

### 1. Input Data File

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


### 2. Annotation File 

The annotation requires a header to specify the name and the nature (categorical or continuous) of each SNP-level anntation. For example,
```
SNP   feature1_d	feature2_c
chr1.51479  0		0.23		
chr1.52058  2		0.45
chr1.52238  1		0.97
```
The first column  (always named "SNP") represents the SNP name. The following columns represent specific annotations. For categorical/discrete annotations, the annotation name should alwasy have a suffix "_d"; whereas for continuous annotations, the name should end with "_c".  In the above example, the annotation "feature1" is a categorical variable with 3 categories (0, 1 and 2), and "feature2" is a continuous variable.There is no restriction on the number of annotations for enrichment analysis. As an example the annotations used for ([Wen et al. 2016, PLoS Genetics](http://www.sciencedirect.com/science/article/pii/S0002929716300957)) can be found [here](http://genome.grid.wayne.edu/centisnps/anno/).

### 3. Gene and SNP Map Files

These two files are typically used in cis-eQTL mappings. The file formats follows MatrixEQTL formats for respective map files. The gene and SNP IDs should be consistent with what used in the input data file.

* The gene map file has the following format 
```
Gene  Chromosome TSS TES
ENSG00000237683  1  139379 139379
ENSG00000237491  1  714162 714162
ENSG00000230021  1  741274 741274
ENSG00000187634  1  860260 860260
...
```
The header is allowed but not required. Note that current implementation only uses transcription start site (TSS) information, the transcription end site (TES) column is currently utilized, therefore, we provide dummy information in above example. 

* The SNP map file has the following format 
```
SNP Chromosome Position
rs1234  1      51485
rs1235  1      52085
...
```
The header is allowed but not required.  



## Output File Format

### 1. Enrichment Estimation File

The enrichment estimation output has the following format

```
Intercept	-8.723	      -8.802	 -8.701	
binding.1	 0.623         0.226      1.021
binding.2	 1.222         0.946      1.499
```
The first column represents the annotation name and its corresponding level (for a categorical variable). The second column is the point estimate (MLE) of the log odds ratio. Columns 3-4 represent the boundries of the corresponding 95% confidence interval.



### 2. QTL Discovery Result

The output from the QTL discovery has the following format
```
    1       ENSG00000164308    3.406e-106    1
    2       ENSG00000166750    4.212e-105    1
    3       ENSG00000198468    1.544e-104    1
    4       ENSG00000174652    6.366e-102    1
    5       ENSG00000197728    1.568e-101    1
    6       ENSG00000233927    1.364e-99    1
    7       ENSG00000237651    5.658e-95    1
    8       ENSG00000006282    8.019e-91    1
    9       ENSG00000112306    5.995e-89    1
   10       ENSG00000166839    9.625e-88    1
```
The output is a ranked list of all tested loci. The first and the second columns represent the rank and the name of a gene/locus, respectively. Column 3 represents the false discovery probability of the corresponding gene (smaller value indicates that the locus likely to harbor QTNs). Finally, column 4 represents a indicator for the hypothesis testing outcome: 1 indicates rejection of the null hypothesis at the pre-defined FDR level.



