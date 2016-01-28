# Multi-SNP Fine-mapping Analysis


## Sample Data Download

A sample data set of gene ENSG00000112799 from the GEUVADIS project is included in the [sample_data](sample_data/) directory.
* [ENSG.dat](sample_data/ENSG.dat)
* [grid](sample_data/grid)
* [ENSG.null.prior](sample_data/.null.prior)
* [ENSG.binding.prior](sample_data/.binding.prior)
 

## Input Data Format

### Genotype-Phenotype Data File (Required)

This data file contains individual-level phenotype and genotype information required by genetic association analysis, and in plain text format. The data file allows information multiple subgroups (e.g. in a meta-analytic setting). Importantly, each data file for multi-SNP fine-mapping analysis should (ideally) be corresponding to a single locus.

The phenotype section contains lines starting with the keyword "pheno". Each line has the following format
```
  pheno  pheno_id subgroup_id pheno_ind_1 pheno_ind_2 ... pheno_ind_n
```  
The pheno\_id field contains a character string that denotes the name of the phenotype. The group\_id field is a character string that uniquely labels a specific group. Note that the group IDs should be consistently used in the data file. The following entries are numerical values of phenotypes for individual 1 to _n_ in the indicated subgroup. Note, because each subgroup can have different number of individuals, the length of each line can differs.

  
The genotype section directly follows the phenotype section. Each line contains the genotype information of a SNP for the samples in a subgroup, i.e.,
```
  geno snp_id  group_id geno_ind_1 geno_ind_2 ... geno_ind_n
```
The leading "geno" is a keyword that indicates the line encodes genotypes. The snp\_id field contains a character string that denotes the ID of a SNP. The additional group_id field indicates the particular subgroup in which genotypes are measured. The remaining entries are genotypes of the target SNP for individual 1 to _n_ in the subgroup coded in dosage format (i.e, 0,1 or 2, or any fractional numbers between [0,2] if the genotype is imputed). 

Note that if the data only contains a single group, the group\_id becomes nuisance. Nevertheless, it is still required.

Missing values in phenotype and genotype data are allowed, and they should be represented by "NA".




### Effect Size Grid File (Required)

The grid file contains information on prior effect size specifications for causal variants. In the fine-mapping analysis, the effect size of each causal variant is treated as a nuisance parameter and integrated out. In our model, the effect size is defined on the scale of signal-noise ratio, and therefore unit-less. For the general meta-analytic setting, two (hyper-)parameters, phi and omega, are required: phi represents the heterogeneity of a genetic effect across multiple subgroups, and omega describes the average effect size. For commonly-used fixed-effect meta-analysis, simply specify phi=0. Alternatively, omega^2/(omega^2 + phi^2) represents the correlation of effects between subgroups, and (omega^2 + phi^2) represents the overall genetic effect which can be linked to heritability (for a SNP with frequency f, the hertiability explained by the single SNP association is roughly 2\*f\*(1-f)\*(omega^2+phi^2)/[1+2\*f\*(1-f)\*(omega^2+phi^2)].


We use a grid of (phi, omega) combinations to describe the long-tail behavior of the genetic effects observed in practice. If the data only contains a single study, only omega needs to be specified and phi should be set as 0. The mathematical details on the effect size prior specification can be found in [Wen and Stephens, 2014](http://projecteuclid.org/euclid.aoas/1396966283) and [Flutre et al,2013](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003486). 

The grid file contains a two-column data matrix: the first column  specifies phi and the second column specifies omega the average effect size (hyper)-parameter (omega). Each row of the grid data matrix provides a unique combination. A sample grid file for fixed-effect meta-analysis is given below.    

```
0.0000  0.1000
0.0000  0.2000
0.0000  0.4000
0.0000  0.8000
0.0000  1.6000
```

The [sample file](sample_data/grid) used in cross-population cis-eQTL analysis allows some effect size heterogeneity for eQTLs across population groups. 


### Prior File (Optional)


#### Automatic Prior Computation by TORUS


### Running Multi-SNP Fine-mapping Analysis

#### Important Command-line Options


### Interpretation of Results











