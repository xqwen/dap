# QTL Discovery using TORUS



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

  


