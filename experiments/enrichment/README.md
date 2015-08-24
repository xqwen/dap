# Simulations for Enrichment Analysis


## 1. Data generation

### Required utilities

* openmp_wrapper []
* sim_data.pl (included)
* sim_pheno.R (included)
* assemble.pl (included)


### Setup

The following files and directories are assumed to be existing by the simulation script:

* sbams_data (empty)
* truth (empty)
* pheno_data (empty)
* geno_data (with existing genotype data, examples in geno_data.tgz with 1500 genes from the GEUVADIS project)
* gene.list (a list of gene names matching genotype data, example included)


### Command

```
perl sim_data.pl lambda_value
```

lambda_value is the true enrichment parameter


### Output

Assembled phenoytpe-genotype data for all genes are formated in sbams format and saved in sbams_data directory, simulated SNP annotations and the simulation truth (e.g. the true QTLs and their effect sizes) are saved in the truth directory.


## 2. Analysis Methods


### (1) Best case analysis 

This analysis assumes the latent association status of all SNPs are known, and simply run a logistic regression to obtain the enrichment parameter estimates.

#### Required utilities

* logit_reg.pl
* logitReg.R

#### Command 

```
perl logit_reg.pl > output
```

#### Output 

1. first column: estimate of the intercept term 
2. second column: point estimate of the enrichment parameter
3. Columns 3-4: 95% CI of the enrichment estimate


### (2) Naive analysis

This analysis carries out a two stage procedure: first classify assocation status of a SNP according to its single SNP association test (FDR level 0.01) and then run logistic regression using the classifcation results.


#### Required utilities

* naive_enrich.pl
* enrich_std.R

#### Command

``` 
perl naive_enrich.pl > output 
```

#### Output

Same as the best case analysis



### (3) Analysis with DAP embedded EM

This analysis runs with DAP embedded EM algorithm. 

 
#### Setup

* Creat an empty directory "dap_rst" before running the analysis


#### Required utlities

* pre-compiled dap executible in the system search path
* em_fmap.pl 
* enrich_m.R

#### Command

For DAP-1 algorithm

```
perl em_fmap.pl -msize 1 

```
 
For adaptive DAP with inclusion threshold 0.05

```
perl em_fmap.pl -it 0.05
```

#### Output

EM iteration information along with the final point estimates and 95% CI of the enrichment parameters are output to screen.


