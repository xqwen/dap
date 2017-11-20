# Exact Bayesian Computation vs. DAP

## 1. Exact Bayesian Computation

This is achieved by running dap with specialized parameters

### command 

``` 
dap -d data_file  -g grid_file -it 0 -st 0  -t nthread
```
The "-it 0" option will force dap to include all SNPs as candidates and the "-st 0" option will force dap to run through all model size partitions, i.e., overide the default stopping rule. We strongly recommend parallel processing using multiple threads with the "-t nthread" option.

### Input files

* data_file:  genotype-phenotype data file in sbams format
* grid_file:  prior effect size required by the Bayes factor calculation


Examples of both files are provided.


### Output

The output is similar to the DAP output, which includes posterior probability for each individual model and PIP for each candidate SNP.


## 2. Simulation 

The simulation procedure takes genotypes from a pre-defined genomic region with 15 SNPs in the GEUVADIS project and simulate quantitative phenotypes for each population. It outputs a genotype-phenotype file in sbams format that is ready for both the exact Bayesian computation and the DAP.


### Required files and directories

* genotype files in directory ``geno_data`` (examples provided in [geno_data.tgz](http://www-personal.umich.edu/~xwen/dap/data/geno_data.tgz))
* pheno_data (empty directory)
* sim_pheno.R 
* sim.pl

### Command

```
perl sim.pl nsig
```

nsig is the number of causal SNPs.

### Output

An output file named ``gene.sbams.dat``


  
