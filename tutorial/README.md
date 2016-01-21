# Tutorial for Integrative Genetic Association Analysis

## Overview of Integrative Genetic Associaiton Analysis


Comparing to the traditional genetic association analysis, which typically attempts to identify association signals between a complex trait and densly genotyped genetic markers (SNPs), the integrative analysis also quantitatively includes genomic annotations of the genetic markers into the association analysis. Our software package aims to address three inter-related analysis goals:

1. Assess the enrichment level of the annotated SNPs in the association signals (Enrichment Analysis)
2. Discover genetic loci that harbor causal variants (QTL Discovery)
3. Perform multi-SNP fine-mapping analysis for the identified loci from 2 (Multi-SNP Fine-mapping)

The first two goals can be achieved by the executable [torus](https://github.com/torus/) and the third aim can be achieved by the executable ``dap``. 



### Types of Applications

We currently support two types of applications: molecular (cis) QTL mapping and tradition single phenotype genome-wide association study (GWAS). In comparison to GWAS,  a distinct feature of molecular QTL mapping is that tens of thousands (or hundreds of thousands) of molecular phenotypes (e.g., gene expression, DNA methylation, chromatin accessibility, histone modifications) are simultaneously measured and analyzed. In addition, the candidate (cis) genomic region for each molecular phenotype is typically not large (usually spanning 1 to 2 Mb), whereas for GWAS the candidate SNPs cover the whole genome. 


### Define Genomic Loci
In molecular QTL mapping, the candidate (cis) locus for each molecular phenotype is naturally defined. For GWAS, we adopt the partitioning algorithm recently proposed by [Berisa and Pickrell, 2015](http://bioinformatics.oxfordjournals.org/content/32/2/283) to segment the whole genome into a set of disjoint loci, which roughly represent independent LD blocks. The partition is population specific: for European population, there are about 1,700 loci and  each locus on average spans 1.6Mb. The detailed information on the partitioning is provided in [here](https://bitbucket.org/nygcresearch/ldetect) by the Pickrell lab.
 

### Supported Data Structure

We currently support genetic association data collected in a single study or in a meta-analytic setting. We are working on extending the software to support applications like multi-tissue eQTL mapping as described in [Flutre et al, 2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003486). 


### Summary Statistics vs. Individual-Level Genetic Data

Currently, the multi-SNP fine-mapping analysis requires individual-level genotype data, and we are actively working to extend the fine-mapping analysis using only summary-level data. 

Both enrichment analysis and QTL discovery require only summary statistics (in the simplest case, z-score or p-value from the single SNP association test). 


 
