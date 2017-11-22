# C++ implementation of Adaptive DAP Algorithm

## Compilation 

The following C/C++ libraries are required for compiling the source code

* GNU [GSL](http://www.gnu.org/software/gsl/) library 
* [OpenMP](http://openmp.org/wp/openmp-compilers/) compiler (many popular compilers are compatible)

Simply run ``make`` to compile the executable named ``dap``


## Command-line Syntax
```
dap -d data_file |  -d_z zvalue_file -d_ld ld_file | -d_est effect_estimate_file -d_ld ld_file -d_n sample_size -d_syy syy [-g grid_file] [-p prior_file] [-msize K] [-it value] [-st value] [-t nthread] [-o output_file] [-l log_file] [--output_all]
```


## Important Tips

Run adaptive DAP algorithm with multi-thread option (```-t nthread```) whenever possible! 

## Input Data Options

DAP now accepts three different types of input data: 
    
* Individual-level genotype and phenotype data in sbams format
* Z-scores for each candidate SNP from single-SNP analysis and an LD matrix
* Estimated effect size (beta-hat) and corresponding standard error for each candidate SNP from single-SNP analysis, an LD matrix, plus two numbers: sample size and Syy (explained below).



