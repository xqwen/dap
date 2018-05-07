# GTEx Whole Blood eQTL Results (v6p)

The file ```gtex_v6p_signal_95_CS.summary``` contains information on eQTL signals from the GTEx whole blood data whose SPIP >= 0.95 (hence a 95% Bayesian credible set can be constructed). The format of the file is 

```
signal_ID   SPIP  size_of_95_credible_set
```

The signal ID is coded by the target gene name and a numerical value. 


The file ```gtex_v6p_signal_95_CS.rst``` provides the detailed information for the 95% credible sets. The file format is 

```
signal_ID SNP_name signal_SPIP snp_PIP
```
Note, only SNPs belong to a credible set is included.

