# QTL Discovery using TORUS



## Data Preparation 

The QTL discovery requires same data format as in the enrichment analysis. It takes summary statistics in forms of [point estimate and standard error](https://github.com/xqwen/dap/tree/master/tutorial/enrichment/qtl/gtex_liver/#input-data-format), [pre-computed Bayes factors](https://github.com/xqwen/dap/tree/master/tutorial/enrichment/qtl/geuvadis/#input-data-format), or simply [z-scores](https://github.com/xqwen/dap/tree/master/tutorial/enrichment/gwas#input-data-format). 

Functional annotations for SNPs and their position information can also be used in the QTL discovery. The data format can be found in above links.


## Running QTL discovery

Running QTL discovery analysis is simply achieved by executing ```torus``` with command line option "-qtl". Here we demonstrate two examples

* QTL discovery in GEUVADIS eQTL data. This analysis uses the annotations of SNP distance
```
torus -d geuv.summary.bf.gz --load_bf -smap geuv.snp.map.gz -gmap geuv.gene.map.gz -qtl > geuv.egene.rst
```

* GWAS analysis of HDL.
```
 torus -d HDL.z_score.gz -annot 1000G.annot.gz --load_zval -qtl > HDL.qtl.rst
````

If no annotations provided, the QTL discovery procedure treats every SNP equally a priori.

### Output from QTL discovery

The output from the QTL discovery has the following format

```  
    1       ENSG00000164308    3.406e-106    1
    2       ENSG00000166750    4.212e-105    1
    3       ENSG00000198468    1.544e-104    1
    4       ENSG00000174652    6.366e-102    1
    5       ENSG00000197728    1.568e-101    1
    6       ENSG00000233927    1.364e-99    1
    7       ENSG00000237651    5.658e-95    1
    8       ENSG00000006282    8.019e-91    1
    9       ENSG00000112306    5.995e-89    1
   10       ENSG00000166839    9.625e-88    1

```
The output is a ranked list of all tested loci. The first and the second columns represent the rank and the name of a gene, respectively. Column 3 represents the false discovery probability of the corresponding gene (smaller value indicates that the locus likely to harbor QTNs). Finally, column 4 represents a indicator for the hypothesis testing outcome: 1 indicates rejection of the null hypothesis at the pre-defined FDR level (5% by default).

  


