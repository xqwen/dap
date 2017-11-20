# C++ implementation of improved Adaptive DAP Algorithm 

## Compilation 

The following C/C++ libraries are required for compiling the source code

* GNU [GSL](http://www.gnu.org/software/gsl/) library 
* [OpenMP](http://openmp.org/wp/openmp-compilers/) compiler (many popular compilers are compatible)

Simply run ``make`` to compile the executable named ``dap-g``


## Major new features

* Automatic signal clustering: the output now automatically identifies potentially multiple independent association signal clusters, and provides information on cluster memberships and the correlation within and between clusters
* Improved efficiency: the adaptive DAP algorithm runs more efficiently by focusing on a relatively smaller set of high-probability models
* Default grid specification: there is a set of default grid values built-in, which is appropriate for analysis for a single genotype-phenotype data set or fixed-effect meta-analysis of multiple data sets.


## Command-line Syntax

```dap-g -d data_file [-g grid_file] [-prior prior_file] [-t nthread] [-o output_file] [-all]```


## Important Tips

Run adaptive DAP algorithm with multi-thread option (``-t nthread``) whenever possible! 


## Command-line Options

### 1. Required input data files

* ``-d data_file``: specify the input genotype-phenotype data file in sbams format ([example](sample_data/sample.sbams.dat) provided)

### 2. Optional input data file

* ``-g grid_file``: specify the prior effect size for Bayes factor calculation ([example](sample_data/sample.fixed_effect.grid) provided)

* ``-prior prior_file``: pre-computed priors for each SNP (example provided), if not specified, the default prior is 1/(number of SNPs), which implies the prior expected association signal is 1. 




### 3. DAP algorithm options

* ``-t nthread``: enable nthread parallel multi-threading computation. This feature is important to achieve the computational efficiency
  

### 4. Output options

* ``-o file_name``: re-direct the output to the specified file. By default, the results output to screen
* ``-all``: output PIP values for all SNPs. By default, only SNPs with high PIPs are shown  


## Input File Format

### Phenotype-Genotype Input 

The input data file for phenotype-genotype data are in "sbams" format. It is a text file designed to represent data from a meta-analytic setting including data from multiple subgroups. 
The file contains a phenotype section and a genotype section.

In the phenotype section, each line contains a list of expression levels of all individuals in a subgroup, i.e.,
```
  pheno  pheno_id group_id exp_ind_1 exp_ind_2 ... exp_ind_n
```
The leading "pheno" is a keyword that indicates the line encodes phenotype measurements. The pheno_id field contains a character string that denotes the name of the phenotype. The group_id field is a character string that uniquely labels a specific subgroup. The following entries are numerical values of expression levels of individual 1 to $n$ for the target gene in the indicated subgroup. Note, because each subgroup my have different number of individuals, the length of each line in this section can differ.


The genotype section directly follows the phenotype section. Each line contains the genotype information of a SNP for samples in a subgroup, i.e.,
```
  geno snp_id  group_id geno_ind_1 geno_ind_2 ... geno_ind_n
```
The leading "geno" is a keyword that indicates the line encodes genotypes. The snp_id field contains a character string that denotes the ID of a SNP. 
The additional group_id field indicates the particular subgroup in which genotypes are measured. The remaining entries are genotypes of the target SNP for individual 1 to n in the subgroup coded in dosage format (i.e, 0,1 or 2, or any fractional number in [0,2] if genotypes are imputed). 

Note that if there is only a single group of data, i.e., not a meta-analysis, the group_id becomes a nuisance. Nevertheless, it is required in the input file.


### Grid File
 
The grid file contains prior specifications for genetic effect sizes. In all cases, the grid file always contains a two-column data matrix: the first column always represents the heterogeneity parameter (phi) and the second column is used to specify the average effect size parameter (omega). Each row of the grid data matrix provides a unique prior model and different rows can be used to describe different prior heterogeneity levels and/or prior average effect sizes.   

In the following sample, the grid assumes five levels of overall prior effects (sqrt(omega^2 + phi^2) = 0.1,0.2,0.4,0.8,1.6) values and a fixed effect across all subgroups. 
```
0.0000  0.1000
0.0000  0.2000
0.0000  0.4000
0.0000  0.8000
0.0000  1.6000
```
If there is only a single group of data, only the overall prior effects, i.e., sqrt(omega^2 + phi^2), affect the analysis. Note that the overall prior effect can be viewed as the prior expected chi^2 value for an associated SNP. 




## Output

The output from the DAP include the following information

* Basic summary information (posterior expected number of QTLs, the full log-likelihood in Bayes factor)

* QTL models with high posterior probabilities

* SNPs with high PIP values

* Signal cluster characterization for high-confidence signals






## Running sample data

``dap-g -d sample.sbams.dat -p sample.prior``

 