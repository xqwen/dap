# Enrichment Analysis of HDL-associated Genetic Variants 

In this example, we study the enrichment of functionally annotated genetic variants associated with high density cholesterol (HDL). In this analysis, we use the z-scores of single-SNP association testing of HDL association ([Pickrell, 2014](http://www.cell.com/ajhg/abstract/S0002-9297(14)00263-8)) and functional annotations derived from ENCODE data by [Gusev et al, 2014](http://www.cell.com/ajhg/abstract/S0002-9297(14)00426-1).  


## Sample Data Download

* z-scores from single-SNP association testing: [HDL.z_score.gz](http://www-personal.umich.edu/~xwen/download/gwas_hdl/HDL.z_score.gz)
* SNP annotation file: [1000G.annot.gz](http://www-personal.umich.edu/~xwen/download/gwas_hdl/1000G.annot.gz)


## Input Data Format

### Z-scores from single-SNP association testing
Z-scores from single-SNP association testing are organized in the following formt
```
   chr1:998395  Loc1  -0.178471
   chr1:1000156  Loc1  -0.169669
   chr1:1001177  Loc1  -0.247359
   chr1:1002932  Loc1  -0.240580
   chr1:1003629  Loc1  -0.169000
   chr1:1004957  Loc1  -1.145393
   chr1:1004980  Loc1  -1.145393
   chr1:1006223  Loc1  -1.174756

```
The first column denotes the SNP IDs, and the second column indicate the LD block of the corresponding SNP. Note, the LD blocks are defined based on the results of [Berisa and Pickrell, 2015](http://bioinformatics.oxfordjournals.org/content/32/2/283). The last column represents the z-scores.



### SNP annotation file

The SNP annotation file contains SNP-level genomic annotations used by TORUS analysis. The annotation file uses a header to specify the number and the nature (categorical or continuous) of the anntations. For example,
```
SNP  annot_d
chr1:226580387 5
chr1:162736463 5
chr1:222359612 0
chr1:157255396 0
chr1:95166832 0
chr1:66857915 0
chr1:63432716 4
chr1:8640831 5
chr1:209894785 5
```
The first column with the header "SNP" represents the SNP name. The following columns represent specific annotations.For categorical/discrete annotations, the header should always have a suffix "_d"; whereas for continuous annotations, the header should ends with "_c". Note that, if a SNP is not annotated (i.e. not appeared) in the annotation file, the default category  0 (i.e., the baseline) is assigned. Nevertheless, we strongly recommend users to annotate all the candidate SNPs.   

In this particular annotation file, the code for the categories represents: 1-coding SNP; 2-utr region; 3-promoter region; 4-DHS region; 5-Intron; 0-baseline/all others.

Finally, all input files should be gzipped. 



## Running Enrichment Analysis

The compiled binary executable ```torus``` is required to run the enrichment analysis. Use the following command to start the enrichment analysis
```
 torus -d HDL.z_score.gz -annot 1000G.annot.gz  -est --load_zval > HDL.enrichment.est
```
In particular,```-est``` instructs TORUS to output the 95% confidence intervals for each estimated enrichment parameter; ```--load_zval``` informs TORUS that the input summary-statistics are z-scores (alternatively, Bayes factors can be pre-computed).


### Output from enrichment analysis

The results for enrichment analysis is directly output to the screen, and can be re-directed to a file (in our example, "gtex_liver.enrichment.est"). The output has the following format
```
Intercept    -11.523        (-11.549, -11.497)
  annot.1      4.684        (2.010,     7.359)
  annot.2      4.585        (0.848,     8.321)
  annot.3      3.940        (2.323,     5.558)
  annot.4      1.688        (0.475,     2.902)
  annot.5      1.550        (0.959,     2.141)
```
The first column represents the annotation name and its corresponding level (for a categorical variable). The second column is the point estimate (MLE) of the log odds ratio. Columns 3-4 represent the 95% confidence interval for the corresponding point estimate. It is a feature of GWAS enrichment analysis that confidence interval can be much larger in comparison to molecular QTL mapping (due to relatively less strong association signals).


