# Simulation Study

This folder contains scripts/code, simulated data and outcomes from our simulation study described in the paper. 

## Data generation

The genotype data ```data_generation``` and the Rscript used to generate the phenotype data are in the directory ```data_generation```.


## Running analysis

+ Run DAP-G with sufficient summary stats: ```analysis_command/batch_dap_sss.cmd````
+ Run DAP-G with z-scores: ```analysis_command/batch_dap_z.cmd```
+ Run FINEMAP with z-scores: ```analysis_command/batch_finemap.cmd```; required data description file ```analysis_command/data```

Note: the input files for FINEMAP and DAP-G with z-scores are identical. For example, ```finemap_data/region1.ld``` is a soft link to ```summary_data/sim.1.LD.dat``` and ```finemap_data/region1.z``` is a soft link to ```summary_data/sim.1.zval.dat```


## Download simulated data and results

+ Simulated individual-level data (in sbams format) and truth [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.sbams_truth.tgz)
+ Summary statistics extracted from individual-level data [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.summary_stats.tgz)
+ Output from FINEMAP [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.finemap_out.tgz)
+ Output from dap-g with sufficient summary statistics [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.dap_out.tgz)
+ Output from dap-g with z-scores [download](http://www-personal.umich.edu/~xwen/dapg_sim/sim_data.dap_z_out.tgz)






