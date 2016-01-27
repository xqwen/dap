#  Cis-eQTL Enrichment Analysis  using GEUVADIS Data 

The GEUVADIS data are multi-population eQTL data, which are originally distributed from [EBI website](http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/). We performed additional pre-processing steps that are documented in [Wen et al, 2015](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005176), and formated the data suitable for our software.




## Sample Data Download

We have prepared the complete input files of GEUVADIS data for user's reference. Due to the size limitation, all the gzipped files are placed on an external server.

* Pre-computed single SNP Bayes factors by sbams_sslr: [geuv.summary.bf.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.summary.bf.gz)

* SNP map: [geuv.snp.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.snp.map.gz)

* Gene map: [geuv.gene.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.gene.map.gz)

* Binding variants annotation file: [geuv.annot.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.annot.gz)



## Input Data Format

The GEUVADIS data contain expression-genotype data from 5 different populations (YRI, CEU, GBR, FIN and TSI), and we perform eQTL analysis across all 5 populations. In this analysis, we consider two types of genomic annotations: SNP distances to TSS and SNP-level transcritpion factor binding annotations from the CENTIPEDE model. 


### Bayes factor input file
We first compute single SNP Bayes factors by software sbams_sslr. The resulting Bayes factor file has the following format:
```
chr20.49475314 ENSG00000000419 -0.087
chr20.49475476 ENSG00000000419 0.077
chr20.49475514 ENSG00000000419 -0.275
```
where the columns 1-3 represent SNP name, gene name and log10 Bayes factor for the corresponding gene-SNP pair.

### SNP and gene map files

To control for SNP distance to TSS, TORUS requires gene TSS information and SNP position files (i.e., gene map and SNP map). 

The gene map  has the MatrixEQTL format,
```
ENSG00000237683  1  139379 139379
ENSG00000237491  1  714162 714162
ENSG00000230021  1  741274 741274
```
where the column 1-3 represent gene name, chromosome and position of TSS. The last column (column 4) is reseved but currently not in use by TORUS (in this example, we just replicate the TSS, but one can certainly replace it with TES information).

The SNP map also has the MatrixEQTL format, i.e., 
```
chr1.51479  1  51479
chr1.52058  1  52058
chr1.52238  1  52238
```
where the column 1-3 represent SNP name, chromosome and SNP position. 

### SNP annotation file

The SNP annotation file contains SNP-level genomic annotations used by TORUS analysis. The annotation file uses a header to specify the number and the nature (categorical or continuous) of the anntations. For example,
```
SNP   binding_d
chr1.51479  0
chr1.52058  2
chr1.52238  1
```
The first column with the header "SNP" represents the SNP name. The following columns represent specific annotations.For categorical/discrete annotations, the header should alwasy have a suffix "_d"; whereas for continuous annotations, the header should ends with "_c". 

In this particular example, category "0" represents the baseline SNPs, category "1" represents SNPs located in DNAseI footprint regions but having little impact on TF binding and category "2" represents SNPs that are computationally predicted to impact TF binding. 


Finally, TORUS expects all input files are gzipped. 


## Running Enrichment Analysis 

The compiled binary executable torus is required to run the enrichment analysis. Use the following command to start the enrichment analysis
```
 torus -d geuv.summary.bf.gz --load_bf -smap geuv.snp.map.gz -gmap geuv.gene.map.gz -annot geuv.annot.gz  -est > geuv.enrichment.est
```
where ``` --load_bf``` specifies the input file is using pre-computed Bayes factors, and ```-est``` instructs TORUS to output the 95% confidence intervals for each estimated enrichment parameter.


## Output and Result

The results for enrichment analysis is directly output to the screen, and can be re-directed to a file (in our example, "geuv.enrichment.est"). The output has the following format
```
binding.1      0.623         0.226      1.021
binding.2      1.222         0.946      1.499
```
The first column represents the annotation name and its corresponding level (for a categorical variable). The second column is the point estimate (MLE) of the log odds ratio. Columns 3-4 represent the 95% confidence interval for the corresponding point estimate.


