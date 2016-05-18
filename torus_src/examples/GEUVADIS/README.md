# Tutorial: eGene Discovery using GEUVADIS Data 

Here, we detail the steps to perform eGene discovery using the real GEUVADIS data. The data set contains multi-population eQTL data and is originally distributed from [this website](http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/analysis_results/). We performed additional pre-processing steps that are documented in [Wen et al, 2015](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005176). 



## Step 1: Prepare input files

The GEUVADIS data contain expression-genotype data from 5 different populations (YRI, CEU, GBR, FIN and TSI), and we perform eQTL analysis across all 5 populations. In this analysis, we consider two types of genomic annotations: SNP distances to TSS and SNP-level transcritpion factor binding annotations. 


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


Finally, TORUS expects all input files are gzipped. 

We have prepared the complete input files of GEUVADIS data for user's reference. Due to the size limitation, all the gzipped files are placed on an external server.

* Pre-computed single SNP Bayes factors by sbams_sslr: [geuv.summary.bf.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.summary.bf.gz)

* SNP map: [geuv.snp.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.snp.map.gz)

* Gene map: [geuv.gene.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.gene.map.gz)

* Binding variants annotation file: [geuv.annot.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.annot.gz)



## Step 2: Enrichment analysis (optional)

The enrichment analysis is embedded in the eGene discovery, it does not need to be run separately. However, if one is interested in obtaining point and uncertainty estimates for the enrichment parameters, use the following command 
```
 torus -d geuv.summary.bf.gz --load_bf -smap geuv.snp.map.gz -gmap geuv.gene.map.gz -annot geuv.annot.gz  -est > geuv.enrichment.est
```
where ``` --load_bf``` specifies the input file is using pre-computed Bayes factors, and ```-est``` instructs TORUS to output the 95% confidence intervals for each estimated enrichment parameter.


### Output from enrichment analysis

The results for enrichment analysis is directly output to the screen, and can be re-directed to a file (in our example, "geuv.enrichment.est"). The output has the following format
```
binding.1      0.623         0.226      1.021
binding.2      1.222         0.946      1.499
```
The first column represents the annotation name and its corresponding level (for a categorical variable). The second column is the point estimate (MLE) of the log odds ratio. Columns 3-4 represent the 95% confidence interval for the corresponding point estimate.



## Step 3: QTL discovery

For QTL discovery accounting for annotations of 1) SNP distance to TSS and 2) Binding variants annotations, issue the following command 
```
torus -d geuv.summary.bf.gz --load_bf -smap geuv.snp.map.gz -gmap geuv.gene.map.gz -annot geuv.annot.gz  -qtl > geuv.egene.rst
```

### Output from QTL discovery

The results for QTL discovery is directly output to the screen, and can be re-directed to a file (in our example, "geuv.egene.est"). The output has the following format

```
1    ENSG00000013573.12    1.616e-224    1
2    ENSG00000163682.11    5.934e-222    1
3    ENSG00000174652.13    1.215e-187    1
4    ENSG00000164308.12    3.318e-182    1
5     ENSG00000233927.4    1.774e-164    1
```
The output is a ranked list of all tested loci. The first and the second columns represent the rank and the name of a gene, respectively. Column 3 represents the false discovery probability of the corresponding gene (smaller value indicates that the locus likely to harbor QTNs). Finally, column 5 represents a indicator for the hypothesis testing outcome: 1 indicates rejection of the null hypothesis at 5% FDR level.

  


