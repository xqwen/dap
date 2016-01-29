#  Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors (DAP)


This repository contains the software implementations for a suite of statistical methods to perform genetic association analysis integrating genomic annotations. These methods are designed to perform rigorous enrichment analysis, QTL discovery and multi-SNP fine-mapping analysis in a highly efficient way. The statistical model and the key algorithm, Deterministic Approximation of Posteriors (DAP), are described in this manuscript. 

The repository includes source code, scripts and necessary data to replicate the results described in the manuscript. A detailed tutorial to guide the users through some specific analysis tasks is also included. 

For questions/comments regarding to the software package, please contact Xiaoquan Wen (xwen at umich dot edu).


## Repository directories

* ``src``: the C/C++ source code for adaptive DAP algorithm (for multi-SNP fine-mapping analysis)

* ``tutorial``: provide detailed descriptions on usage of the software package (including data formating, results interpretation etc.)

* ``utility``: utility scripts for results interpretation etc.

* ``experiments``: this directory contains necessary scripts/code and data for evaluating the DAP, e.g., performance comparison with the MCMC, exact Bayesian calculation etc. 


In addition, the source code implementing enrichment analysis and QTL discovery can be found in a separate repo [torus] (https://github.com/xqwen/torus)


## [Tutorial](tutorial/)

We provide a detailed tutorial to demonstrate the usage of software in enrichment analysis, QTL discovery and multi-SNP fine-mapping.


## Contributors

* Xiaoquan Wen (University of Michigan)
* Roger Pique-Regi (Wayne State University)
* Yeji Lee (University of Michigan)

