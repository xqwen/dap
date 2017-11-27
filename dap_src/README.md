# C++ implementation of Adaptive DAP Algorithm

## Compilation 

The following C/C++ libraries are required for compiling the source code

* GNU [GSL](http://www.gnu.org/software/gsl/) library 
* [OpenMP](http://openmp.org/wp/openmp-compilers/) compiler (many popular compilers are compatible)

Simply run ``make`` to compile the executable named ``dap``


## Command-line Syntax
```
dap -d data_file |  -d_z zvalue_file -d_ld ld_file | -d_est effect_estimate_file -d_ld ld_file -d_n sample_size -d_syy syy [-g grid_file] [-p prior_file] [-msize K] [-it value] [-st value] [-t nthread] [-o output_file] [-l log_file] [--output_all]
```


## Important Tips

Run adaptive DAP algorithm with multi-thread option (```-t nthread```) whenever possible! 

## Input Data Options

DAP now accepts three different types of input data: 
    
* Individual-level genotype and phenotype data in sbams format. For this option, use ```-d sbams_data_file``` to speify the single input file.
* Z-scores for each candidate SNP from single-SNP analysis and an LD matrix. Two input files are required for this option, use ```-d_z zval_file  -d_ld LD_file``` to specify the z-score and LD (genotype correlation) files, respectively.
* Estimated effect size (beta-hat) and corresponding standard error for each candidate SNP from single-SNP analysis, an LD matrix, plus two numbers: sample size and total sum of squares of the outcome variable (denoted by Syy). The two input files containing single SNP association estimates and LD information are specified by ```-d_est estimate_file -d_ld LD_file```, respectively; use ```-d_n sample_size -d_syy Syy``` to provide numeric values for the sample size and Syy.  

## Input Data Format

### Individual-level Data 

The input data file for phenotype-genotype data are in "sbams" format. It is a text file originally designed to represent data from a meta-analytic setting including data from multi
ple subgroups. Note that, the *current* implementation of the DAP algorithm only supports a single subgroup (the multi-group feature is maintained in [version 1](../version1/)). We will add the multi-group feature back in the near feature.

The file contains a phenotype section, a controlled covariate section and a genotype section.

In the phenotype section, each line contains a list of measured traits of all individuals in a subgroup, i.e.,
```
pheno  pheno_id group_id exp_ind_1 exp_ind_2 ... exp_ind_n
```
The leading "pheno" is a keyword that indicates the line encodes phenotype measurements. The pheno_id field contains a character string that denotes the name of the phenotype. The group_id field is a character string that uniquely labels a specific subgroup. The following entries are numerical values of expression levels of individual 1 to n for the target gene in the indicated subgroup. Note, because each subgroup my have different number of individuals, the length of each line in this section can differ

```
geno snp_id  group_id geno_ind_1 geno_ind_2 ... geno_ind_n
```
The leading "geno" is a keyword that indicates the line encodes genotypes. The snp_id field contains a character string that denotes the ID of a SNP. The additional group_id field indicates the particular subgroup in which genotypes are measured. The remaining entries are genotypes of the target SNP for individual 1 t n in the subgroup coded in dosage format (i.e, 0,1 or 2, or any fractional number in [0,2] if genotypes are imputed).

Similarly, the controlled covariates are coded in the similar format lead by the keyword "controlled", i.e.,
```
controlled variable_name group_id value_ind_1 value_ind_2 ... value_ind_n
```




