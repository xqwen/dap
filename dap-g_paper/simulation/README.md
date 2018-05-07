# Simulation Study

This folder contains scripts/code, simulated data and outcomes from our simulation study described in the paper. 

## Data generation

The genotype data ```sim.r.geno``` and the Rscript used to generate the phenotype data are in the directory ```data_generation```.

We set up a simulation scenario mimickingcis-eQTL mapping in a practical setting. In particular, we use the real genotype data from 343 European individuals from the GUEVADIS project. 
We artificially construct a genomic region of 1,001 SNPs. The region is divided into 91 LD blocks, and each block contains 11 SNPs.  All LD blocks are selected from chromosome 1, and the consecutive blocks are at least 1Mb apart. 
With this construction scheme, the LD only presents within each block, and the SNP genotypes are mostly uncorrelated across blocks



## Running analysis

+ Run DAP-G with sufficient summary stats: ```analysis_command/batch_dap_sss.cmd```
+ Run DAP-G with z-scores: ```analysis_command/batch_dap_z.cmd```
+ Run FINEMAP with z-scores: ```analysis_command/batch_finemap.cmd```; required data description file ```analysis_command/data```

Note: the input files for FINEMAP and DAP-G with z-scores are identical. For example, ```finemap_data/region1.ld``` is a soft link to ```summary_data/sim.1.LD.dat``` and ```finemap_data/region1.z``` is a soft link to ```summary_data/sim.1.zval.dat```


## Download simulated data and results

+ Simulated individual-level data (in sbams format) and truth [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.sbams_truth.tgz)
+ Summary statistics extracted from individual-level data [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.summary_stats.tgz)
+ Output from FINEMAP [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.finemap_out.tgz)
+ Output from dap-g with sufficient summary statistics [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.dap_out.tgz)
+ Output from dap-g with z-scores [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.dap_z_out.tgz)

## Evaluation of simulation results

The following sub-directories contain code and results for the following analyses

+ Construction of ROC curves: ```ROC```
+ Calibration of SNP-level PIPs: ```calibration```
+ Signal-level FDR control: ```fdr```





