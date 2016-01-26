# Cis-eQTL Enrichment Analysis using GTEx Liver Data 

The data set contains eQTL data in liver which are originally downloaded from [GTEx portal](http://gtexportal.org/home/). In particular, the summary-level data for each gene-SNP pair are directly output from the software package [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/), which is commonly used for cis-eQTL analysis. In this example, we show that the enrichment analysis can take the input and output files from the popular eQTL analysis software ``MatrixEQTL`` and perform the desired enrichment analysis. 


## Sample Data Download

* single-SNP eQTL summary statistics output from ```MatrixEQTL```: gtex.liver.summary.gz (This file is about 3Gb in size and is too large for our server, please download it directly from [GTEx Portal](http://gtexportal.org/home/))
* SNP map: [gtex.snp.map.gz](http://www-personal.umich.edu/~xwen/download/gtex_liver/gtex.snp.map.gz)
* Gene map: [gtex.gene.map.gz](http://www-personal.umich.edu/~xwen/download/gtex_liver/gtex.gene.map.gz)
* Binding variants annotation file from the CENTIPEDE model: [gtex.centipede.annot.gz](http://www-personal.umich.edu/~xwen/download/gtex_liver/gtex.centipede.annot.gz) 


## Input Data Format

### MatrixEQTL output of summary-level statistics in single-SNP testing
The output from MatrixEQTL contains the association information of each gene-SNP pair with the following format
```
SNP	gene	beta	 t-stat	p-value
1_30923_G_T_b37	ENSG00000227232.4	0.00399184253633977	0.0232610959574911	0.981503781647425
1_51479_T_A_b37	ENSG00000227232.4	-0.0264207753473278	-0.101592859539226	0.919350965625195
1_55299_C_T_b37	ENSG00000227232.4	-0.224656599916996	-0.977243563498281	0.331590265844854

```
The columns 1-2 represent SNP and gene names, respectively. Column 3 represents the least-square estimate (LSE, in this case also MLE) of the genetic effect for each gene-SNP pair. Columns 4 and 5 represent the corresponding t-statistic and p-values. 

This output file from MatrixEQTL can be directly taken as input for ``torus``.


### SNP and gene map files

To control for SNP distance to TSS, TORUS requires gene TSS information and SNP position files (i.e., gene map and SNP map). 

The gene map also has the MatrixEQTL format,
```
Id	Chr	TSS	TSS
ENSG00000223972.4	1	11869	11869
ENSG00000243485.2	1	29554	29554
ENSG00000227232.4	1	29806	29806
ENSG00000237613.2	1	36081	36081
```
where the column 1-3 represent gene name, chromosome and position of TSS. The last column (column 4) is reseved but currently not in use by TORUS (in this example, we just replicate the TSS, but one can certainly replace it with TES information).

The SNP map also has the MatrixEQTL format, i.e., 
```
ID	CHROM	POS
1_30923_G_T_b37	1	30923
1_51479_T_A_b37	1	51479
1_52058_G_C_b37	1	52058
1_52238_T_G_b37	1	52238
1_54490_G_A_b37	1	54490
```
where the column 1-3 represent SNP name, chromosome and SNP position. 

### SNP annotation file

The SNP annotation file contains SNP-level genomic annotations used by TORUS analysis. The annotation file uses a header to specify the number and the nature (categorical or continuous) of the anntations. For example,
```
SNP binding_d
1_66162_A_T_b37 1
1_66176_T_A_b37 1
1_66219_A_T_b37 1
1_66331_A_C_b37 1
1_66442_T_A_b37 1
1_66457_T_A_b37 1
1_66507_T_A_b37 1
1_120983_C_T_b37 2
1_121009_C_T_b37 1
```
The first column with the header "SNP" represents the SNP name. The following columns represent specific annotations.For categorical/discrete annotations, the header should always have a suffix "_d"; whereas for continuous annotations, the header should ends with "_c". Note that, if a SNP is not annotated (i.e. not appeared) in the annotation file, the default category  0 (i.e., the baseline) is assigned. Nevertheless, we strongly recommend users to annotate all the candidate SNPs.   

Finally, TORUS expects all input files are gzipped. 



## Running Enrichment Analysis

The compiled binary executable ```torus``` is required to run the enrichment analysis. Use the following command to start the enrichment analysis
```
 torus -d  gtex.liver.summary.gz  -smap gtex.snp.map.gz -gmap gtex.gene.map.gz -annot gtex.centipede.annot.gz  -est > gtex_liver.enrichment.est
```
In particular,```-est``` instructs TORUS to output the 95% confidence intervals for each estimated enrichment parameter.


### Output from enrichment analysis

The results for enrichment analysis is directly output to the screen, and can be re-directed to a file (in our example, "gtex_liver.enrichment.est"). The output has the following format
```
binding.1      0.623         0.226      1.021
binding.2      1.222         0.946      1.499
```
The first column represents the annotation name and its corresponding level (for a categorical variable). The second column is the point estimate (MLE) of the log odds ratio. Columns 3-4 represent the 95% confidence interval for the corresponding point estimate.


