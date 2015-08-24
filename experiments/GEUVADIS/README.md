# GEUVADIS Analysis


## Input files
Here we provide necessary data to carry out the enrichment analysis using the GEUVADIS data with the DAP1 embedded EM algorithm. Due to the size limitation, all the gzipped files are placed on an external server.

* Pre-computed single SNP Bayes factors by sbams_sslr : [geuv.summary.bf.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.summary.bf.gz]

* SNP position file: [geuv.snp.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.snp.map.gz)

* Gene TSS position file: [geuv.gene.map.gz](http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.gene.map.gz)

* Binding variants annotation file: [http://www-personal.umich.edu/~xwen/dap/data/geuv/geuv.annot.gz]


## Running analysis


```
em_dap1 -d geuv.summary.bf.gz --load_bf -smap geuv.snp.map.gz -gmap geuv.gene.map.gz -annot geuv.annot.gz
```


