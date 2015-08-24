# Simulations for Enrichment Analysis


## Data generation

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

Assembled phenoytpe-genotype data for all genes are formated in sbams format and saved in sbams_data directory, simulated SNP annotations and the simulation truth (e.g. the true QTLs and their effect sizes) are saved in the truth directory





