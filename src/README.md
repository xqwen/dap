# C++ implementation of Adaptive DAP Algorithm

## Compilation 

The following C/C++ libraries are required for compiling the source code

* GNU GSL library 
* OpenMP library  

Simply run ``make`` to compile the executable named ``dap``


## Important command line options

### 1. Required input data files

* ``-d data_file``: specify the input genotype-phenotype data file in sbams format (example provided)
* ``-g grid_file``: specify the prior effect size for Bayes factor calculation (example provided)


### 2. Optional input data file

* ``-prior prior_file``: pre-computed priors for each SNP (example provided)




### 3. DAP algorithm options

* ``-msize K``: this option option will enable DAP-K algorithm instead of the adaptive DAP.   
* ``-it value``: the including threshold for screening high priority SNPs. Valid value should be in [0,1). When it is set to 0, all SNPs will be included as high priority  (default value 0.01)
* ``-st value``: additional threshold to determine the model size partitions for extensive exploration. Valid value should be greater or equal to 0. Small value will force the DAP to explore higher size partitions and ensure numerical accuracy.  There are two extremes: 1) if it is set to 0, DAP will igonre the built-in stopping rule by combinatorial approximation and explore all the partitions; 2) if it is set to a large value (e.g. 100), only the built-in stopping rule takes effect. By default, it is set to 0.01.
* ``-t nthread``: enable nthread parallel multi-threading computation. This feature is important to achieve the computational efficiency
  

### 4. Output options

* ``-o file_name``: re-direct the output to the specified file. By default, the results output to screen
* ``-all``: output PIP values for all SNPs. By default, only SNPs with high PIPs are shown  


## Output

The output from the DAP include the following information

* Basic summary information (posterior expected number of QTLs, the full log-likelihood in Bayes factor)
* QTL models with high posterior probabilities
* SNPs with high PIP values


## Running sample data

``dap -d sample.sbams.dat -g sample.fixed_effect.grid -p sample.prior``

 