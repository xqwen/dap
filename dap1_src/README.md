# Fine-mapping with DAP-1 using Summary-level Statistics


```dap1``` implements the DAP-1 algorithm for fine-mapping analysis. The main advantage of the algorithm is that it utilizes only summary-level information (e.g., pre-computed Bayes factors or z-scores) from single-SNP association testing.

Warning: this particular implementation is under active development, use the current code with cautions.

## Compilation 


Simply run ``make`` to compile the executable named ``dap_ss``.


## Command-line Syntax

```dap1 -d data_file [-prior prior_file] > output_file```


## Command-line Options

### 1. Required input data files

* ``-d data_file``: specify the summary-level statistics for a genomic locus of interest



### 2. Optional input data file

* ``-prior prior_file``: pre-computed priors for each SNP (example provided)



## Input File Format


Currently, ```dap_ss``` supports following input formats

* Summary-level Z-statistics from single SNP association tests, e.g.,
```
SNP	locus	z-val
rs1234	geneA	1.28
...
```
A header is allowed (but not required) in the input file. 



## Output

The output has the following format
```
SNP   PIP    log10BF
  chr6:34592090    6.355e-01     6.901
  chr6:34603691    5.646e-02     6.704
  chr6:34616322    3.380e-02     6.482
  chr6:34602685    1.548e-02     5.990
  chr6:34800435    1.173e-02     6.286
  chr6:34557246    9.589e-03     5.773
  chr6:34552797    8.833e-03     6.242
```

## Running sample data

``dap1 -d sample.zval.dat  -p sample.prior``

 