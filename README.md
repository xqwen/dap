#  Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors (DAP)

The current implementation is DAP-G!

This repository contains the software implementations for a suite of statistical methods to perform genetic association analysis integrating genomic annotations. These methods are designed to perform rigorous enrichment analysis, QTL discovery and multi-SNP fine-mapping analysis in a highly efficient way. The statistical model and the key algorithm, Deterministic Approximation of Posteriors (DAP), are described in this [manuscript](http://biorxiv.org/content/early/2016/03/26/026450). 

The repository includes source code, scripts and necessary data to replicate the results described in the manuscript. A detailed tutorial to guide the users through some specific analysis tasks is also included. 

For questions/comments regarding to the software package, please contact Xiaoquan (William) Wen (xwen at umich dot edu).



## License

Software distributed under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See [LICENSE](http://www.gnu.org/licenses/gpl-3.0.en.html) for more details.





## Repository directories

* ``dap_src``: C/C++ source code of the adaptive DAP algorithm with new and improved features (now working with summary-level statistics)


* ``torus_src``: C/C++ source code of the EM-DAP1 algorithm (for enrichment analysis and QTL discovery)


* ``utility``: utility scripts for results interpretation, file format conversion etc.

* ``version 1``: legacy code of the DAP implementation from version 1 

## [Tutorial](https://github.com/xqwen/dap/wiki)

We provide a detailed tutorial wiki to demonstrate the usage of software in enrichment analysis, QTL discovery and multi-SNP fine-mapping. Sample data are provided when possible. 






## Contributors

* Xiaoquan Wen (University of Michigan)
* Roger Pique-Regi (Wayne State University)
* Yeji Lee (University of Michigan)

## Citation

* Wen, X., Lee, Y., Luca, F., Pique-Regi, R. [Efficient Integrative Multi-SNP Association Analysis using Deterministic Approximation of Posteriors.](http://biorxiv.org/content/early/2016/03/26/026450)  The American Journal of Human Genetics, 98(6), 1114-1129
* Lee, Y, Luca, F, Pique-Regi, R,Wen, X. [Bayesian Multi-SNP Genetic Association Analysis: Control of FDR and Use of Summary Statistics](https://www.biorxiv.org/content/early/2018/05/08/316471) biRxiv:316471

