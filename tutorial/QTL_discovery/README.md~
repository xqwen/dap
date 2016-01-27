# Tutorial: eGene Discovery using GTEx Liver Data 

Here, we detail the steps to perform eGene discovery using the GTEx Liver data. The data set contains multi-population eQTL data and is originally distributed from [GTEx portal](http://gtexportal.org/home/). In particular, the summary-level data for each gene-SNP pair are provided in the output format of software package [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/).
The implication here is that if the user has already run MatrixEQTL for single SNP analysis, the results can be directly ported into the ```torus``` for QTL discovery. 




## Step 1: Prepare input files

### MatrixEQTL input file for summary-level statistics
The output from the MatrixEQTL contains the association information of each gene-SNP pair with the following format
```
SNP	gene	beta	 t-stat	p-value
1_30923_G_T_b37	ENSG00000227232.4	0.00399184253633977	0.0232610959574911	0.981503781647425
1_51479_T_A_b37	ENSG00000227232.4	-0.0264207753473278	-0.101592859539226	0.919350965625195
1_55299_C_T_b37	ENSG00000227232.4	-0.224656599916996	-0.977243563498281	0.331590265844854

```
The columns 1-2 represent SNP and gene names, respectively. Column 3 represents the least-square estimate (LSE, in this case also MLE) of the genetic effect for each gene-SNP pair. Columns 4 and 5 represent the corresponding t-statistic and p-values. 

### SNP and gene map files

To control for SNP distance to TSS, TORUS requires gene TSS information and SNP position files (i.e., gene map and SNP map). 

The gene map also has the MatrixEQTL format,
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

The required input files of GTEx liver eQTL data can be downloaded from the [GTEx Portal](http://www.gtexportal.org/home/).


## Step 2: Enrichment analysis (optional)

The enrichment analysis is embedded in the eGene discovery, it does not need to be run separately. However, if one is interested in obtaining point and uncertainty estimates for the enrichment parameters, use the following command 
```
 torus -d  Liver_Analysis.cis.eqtl.gz  -smap snp.map.gz -gmap gene.map.gz -annot gtex.annot.gz  -est > gtex_liver.enrichment.est
```
In particular,```-est``` instructs TORUS to output the 95% confidence intervals for each estimated enrichment parameter.


### Output from enrichment analysis

The results for enrichment analysis is directly output to the screen, and can be re-directed to a file (in our example, "gtex_liver.enrichment.est"). The output has the following format
```
binding.1      0.623         0.226      1.021
binding.2      1.222         0.946      1.499
```
The first column represents the annotation name and its corresponding level (for a categorical variable). The second column is the point estimate (MLE) of the log odds ratio. Columns 3-4 represent the 95% confidence interval for the corresponding point estimate.



## Step 3: QTL discovery

For QTL discovery accounting for annotations of 1) SNP distance to TSS and 2) Binding variants annotations, issue the following command 
```
torus -d  Liver_Analysis.cis.eqtl.gz  -smap snp.map.gz -gmap gene.map.gz -annot gtex.annot.gz   -qtl > gtex_liver.egene.rst
```

### Output from QTL discovery

The results for QTL discovery is directly output to the screen, and can be re-directed to a file (in our example, "gtex_liver.egene.est"). The output has the following format

```
1    ENSG00000013573.12    1.616e-224    1
2    ENSG00000163682.11    5.934e-222    1
3    ENSG00000174652.13    1.215e-187    1
4    ENSG00000164308.12    3.318e-182    1
5     ENSG00000233927.4    1.774e-164    1
```
The output is a ranked list of all tested loci. The first and the second columns represent the rank and the name of a gene, respectively. Column 3 represents the false discovery probability of the corresponding gene (smaller value indicates that the locus likely to harbor QTNs). Finally, column 5 represents a indicator for the hypothesis testing outcome: 1 indicates rejection of the null hypothesis at 5% FDR level.

  


