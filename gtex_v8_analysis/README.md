# Cis-eQTL Fine-mapping of GTEx v8 Data

This directory documents the procedures to fine-map cis-eQTLs in GTEx v8 data. The relevant scripts are also included. For privacy reasons, we are not allowed to distribute individual-level geneotype data. The access to these data must go through dbGAP. 


## Converting to sbams format

1. Run process script
```
perl process.pl -e tissue_expression -g master_genotype_vcf -c tissue_covariate -t tissue_name
```
The ``tissue_expression``, ``master_genotype_vcf``, and ``tissue_covariate`` files are distributed by GTEx. 

In the end, depending on the tissue name, a directory named as ``tissue_name`` and a single batch command file ``tissue_name.assemble.cmd`` are  generated in the current working directory. 

2. Execute batch assemble commands to obtain SBMAS files in ``tissue_name`` directory. This will require script ``assemble.pl`` which is distributed here. 


## Estimate priors for fine-mapping

We directly take the GTEx distributed single-SNP testing output to estimate the fine-mapping priors. In particular, this computation takes into account of SNP distance to transcription start site (TSS). The command to run is
```
 torus -d fastqtl_single_snp_output --fastqtl -dump_prior tissue_name/prior
```
Note that ``tissue_name/prior`` is a directory and should be created manually before running [``torus``](https://github.com/xqwen/torus/). 

If you are not sure about input data format, make sure the following columns are presented in the correct order

    + column 1: gene name
    + column 2: SNP name
    + column 3: distance to TSS
    + column 4: p-value
    + column 5: beta-hat
    + column 6: standard error of beta-hat


For each gene, a prior file will be generated in the ``tissue_name/prior`` directory.


## Running DAP

The command we use for GTEx analysis is 
```
dap-g -d gene_sbams_file -p gene_prior_file -ld_control 0.5 --all -t 4 > output_file
```


``-d`` and ``-p`` point to the sbams and prior files. Option ``--all`` forces outputting information of all SNPs (not just noteworthy ones).  Note that ``--all`` has  now become a default option in DAP. ``-t 4`` indicates the run will use 4 parallel threads if available. ``-ld_control 0.5`` specifies the lowest LD threshold (R^2) to admit a SNP into a signal cluster. 











