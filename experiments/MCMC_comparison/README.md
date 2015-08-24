# DAP vs. MCMC 

The directory contains simulated data and code to compare the performance of MCMC and DAP in mapping multiple QTLs. 

## Simulated data set


### Data files

The file `` simulated_data.smbas.tgz`` contains 1,500 simulated data sets using the genotypes from the GEUVADIS (1000Genome) project. The data sets are all formulated in the sbams format and ready to be processed by the DAP and the MCMC algorithm implemented in the C++ program ``sbams_sslr``.

Each data set contains 2,500 SNPs from the cis-regions of 100 random genes (i.e., 25 neighboring SNPs per region). As a result, the assembled genomic region consists of 100 relatively independent LD blocks with modest to high LD within each block. In each simulation, we randomly assign 1 to 4 causal QTLs and simulate a quantitative trait. The goal of the analysis is to identify the LD blocks that harbor the true QTLs. The truth is recorded in the file ``simulation_trhth.tgz``.




## Analysis methods


### DAP 

``dap -d sim_data.sbams.dat -g grid_file -it thresh -t nthread > dap.out``

* ``sim_data.sbams.dat`` is a simulated data set in sbams format
* ``grid_file`` is necessary for Bayes factro computation and an example is provided


### MCMC

The MCMC algorithm is implemented in the software [sbams_sslr](https://github.com/xqwen/fmeqtl/tree/master/src/sbams_sslr). The command to perform multi-SNP QTL mapping is 
```
sbams_sslr -d sim_data.sbams.dat -g grid_file -b 25000 -r 50000 -meta > mcmc.out
```
The "-b" and "-r" options specify the burnin and repeat lengths for the MCMC run, and the "-meta" option will allow MCMC to perform meta-analysis across population groups.     

