# Tutorial for Integrative Genetic Association Analysis


* [Overview](#overview-of-integrative-genetic-association-analysis)
  * [Types of Applications](#types-of-applications)
  * [Defining Genomic Loci](#defining-genomic-loci)
  * [Supported Data Structure](#supported-data-structure)
  * [Summary Statistics vs. Individual-Level Data](#summary-statistics-vs-individual-level-data)

* [Case Studies](#case-studies)
  * [Enrichment Analysis](#enrichment-analysis)
    * [Enrichment Analysis for cis-eQTL](#enrichment-analysis-in-cis-eqtl-mapping)
    * [Enrichment Analysis for GWAS](#enrichment-analysis-in-gwas)
  * [QTL Discovery](#qtl-discovery)
    * [QTL Discovery in cis-eQTL Mapping](#qtl-discovery-in-cis-eqtl-mapping)
  * [Multi-SNP Fine-mapping](#multi-snp-fine-mapping)


## Overview of Integrative Genetic Association Analysis


Comparing to the traditional genetic association analysis, which typically attempts to identify association signals between a complex trait and densely genotyped genetic markers (SNPs), the integrative analysis also quantitatively includes genomic annotations of the genetic markers into the association analysis. Our software package aims to address three inter-related analysis goals:

1. Assess the enrichment level of the annotated SNPs in the association signals (Enrichment Analysis)
2. Discover genetic loci that harbor causal variants (QTL Discovery)
3. Perform multi-SNP fine-mapping analysis for the identified loci from 2 (Multi-SNP Fine-mapping)

The first two goals can be achieved by the executable "[torus](https://github.com/xqwen/torus/)" and the third aim can be achieved by the executable "dap". 



### Types of Applications

We currently support two types of applications: molecular (cis) QTL mapping and tradition single phenotype genome-wide association study (GWAS). In comparison to GWAS,  a distinct feature of molecular QTL mapping is that tens of thousands (or hundreds of thousands) of molecular phenotypes (e.g., gene expression, DNA methylation, chromatin accessibility, histone modifications) are simultaneously measured and analyzed. In addition, the candidate (cis) genomic region for each molecular phenotype is typically not large (usually spanning 1 to 2 Mb), whereas for GWAS, the candidate SNPs cover the whole genome. 


### Defining Genomic Loci
In molecular QTL mapping, the candidate (cis) locus for each molecular phenotype is naturally defined. For GWAS, we adopt the partitioning algorithm recently proposed by [Berisa and Pickrell, 2015](http://bioinformatics.oxfordjournals.org/content/32/2/283) to segment the whole genome into a set of disjoint loci, which roughly represent independent LD blocks. The partition is population specific: for European population, there are about 1,700 loci and  each locus on average spans 1.6Mb. The detailed information on the partitioning is provided in [here](https://bitbucket.org/nygcresearch/ldetect) by the Pickrell lab.
 

### Supported Data Structure

We currently support genetic association data collected in a single study or in a meta-analytic setting. We are active working on extending the software to support applications like multi-tissue eQTL mapping as described in [Flutre et al, 2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003486). 


### Summary Statistics vs. Individual-Level Data

In the current implementation, the multi-SNP fine-mapping analysis requires individual-level genotype data, and we are actively working to extend the fine-mapping analysis using only summary-level data.

Both enrichment analysis and QTL discovery require only summary statistics (in the simplest case, z-score or p-value from the single SNP association test). 


## Case Studies


### Enrichment Analysis 

#### [Enrichment Analysis in cis-eQTL Mapping](https://github.com/xqwen/dap/tree/master/tutorial/enrichment/qtl/)

This case study provides an example to perform enrichment analysis in molecular QTL mapping. The question we are asking in this case is: is SNP predicted to disrupt transcription factor (TF) binding enriched in cis-eQTLs? We perform the analysis using two eQTL data sets

* [a single-tissue eQTL study using GTEx liver data](enrichment/qtl/gtex_liver/)
* [a cross-population eQTL study using GEUVADIS data](enrichment/qtl/geuvadis/)

In both examples, we use the TF binding annotations from the CENTIPEDE model and account for the genomic position of each candidate SNP with respect to the transcription start site (TSS) of the corresponding target gene.

In our demonstration, we use two types of input format to run the enrichment analysis: the single-tissue analysis of GTEx liver tissue uses the summary-level output from the software package [MatrixEQTL](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/); for the GEUVADIS data, we pre-process the data using a Bayesian single-SNP meta-analysis method and use the resulting Bayes factors as the input for the enrichment analysis. 


#### [Enrichment Analysis in GWAS](https://github.com/xqwen/dap/tree/master/tutorial/enrichment/gwas) 

In this example, we demonstrate the enrichment analysis in GWAS of a complex trait (HDL) using the genomic annotations used in [Gusev et al, 2014](http://www.cell.com/ajhg/abstract/S0002-9297(14)00426-1) (downloaded from [here](https://data.broadinstitute.org/alkesgroup/ANNOTATIONS/)). In this analysis, we use the summary-level single SNP association testing z-scores originally used by [Pickrell, 2014](http://www.cell.com/ajhg/abstract/S0002-9297(14)00263-8). Follow the [link](enrichment/gwas/) for details.

### QTL Discovery 



### Multi-SNP Fine-mapping
