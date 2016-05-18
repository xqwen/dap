#  Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors (DAP)


This repository contains the software implementations for a suite of statistical methods to perform genetic association analysis integrating genomic annotations. These methods are designed to perform rigorous enrichment analysis, QTL discovery and multi-SNP fine-mapping analysis in a highly efficient way. The statistical model and the key algorithm, Deterministic Approximation of Posteriors (DAP), are described in this [manuscript](http://biorxiv.org/content/early/2016/03/26/026450). 

The repository includes source code, scripts and necessary data to replicate the results described in the manuscript. A detailed tutorial to guide the users through some specific analysis tasks is also included. 

For questions/comments regarding to the software package, please contact Xiaoquan Wen (xwen at umich dot edu).




## [Tutorial](https://github.com/xqwen/dap/wiki)

We provide a detailed tutorial wiki to demonstrate the usage of software in enrichment analysis, QTL discovery and multi-SNP fine-mapping. Sample data are provided when possible. 



## Repository directories

* ``dap_src``: C/C++ source code of the adaptive DAP and DAP-K algorithm (for multi-SNP fine-mapping analysis)

* ``torus_src``: C/C++ source code of the EM-DAP1 algorithm (for enrichment analysis and QTL discovery)

* ``sample_data``: sample data for multi-SNP fine-mapping

* ``utility``: utility scripts for results interpretation, file format conversion etc.

* ``experiments``: this directory contains necessary scripts/code and data for evaluating the DAP in the manuscript, e.g., performance comparison with the MCMC, exact Bayesian calculation etc. Note, some of the code has become obsolete, i.e., better code has been implemented in either ```torus``` or ```dap```. This particular directory is setup mostly for reproduciblility purpose.   



## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.


## Contributors

* Xiaoquan Wen (University of Michigan)
* Roger Pique-Regi (Wayne State University)
* Yeji Lee (University of Michigan)

## Citation

* Wen, X., Lee, Y., Luca, F., Pique-Regi, R. [Efficient Integrative Multi-SNP Association Analysis using Deterministic Approximation of Posteriors.](http://biorxiv.org/content/early/2016/03/26/026450)  *bioRxiv* 026450.

