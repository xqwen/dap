### An R package for Bayesian model selection via adaptive deterministic approximation of posteriors

This package is designed for structured Bayesian model selection via DAP algorithm. Applications include genetic association analysis integrating genomic annotations. These methods are designed to perform rigorous enrichment analysis, QTL discovery and multi-SNP fine-mapping analysis in a highly efficient way.

1. Installation
---------------

Pre-installation of [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) is required, and [OpenMP](https://www.openmp.org) compatible compiler is recommended.

`dap` has not been released on CRAN yet. It can be installed as follows in an R session:

``` r
devtools::install_github("https://github.com/yukt/dap.git", ref="rpackage", subdir = "dap_src/R-package")
```

Once installed, the package can be loaded in a given R session using:

``` r
library(dap)
```

2. Analysis Examples
--------------------

To illustrate how the package works, we create a toy dataset `df`. It consists of a response variable `y` and 1000 candidate predictors, in which only `a1` and `b1` are causal. `a2` and `a3` are highly correlated with `a1`, `b2` is highly correlated with `b1`, and `x1`-`x995` are unrelated to `y`.

``` r
set.seed(0)
n = 100
p = 1000

a1 = rnorm(n)
a2 = 10*a1+rnorm(n)
a3 = 2*a1+9*a2+rnorm(n)

b1 = rnorm(n)
b2 = 8*b1+rnorm(n)

x = matrix(rnorm(n*(p-5)), nrow=n)

df = data.frame(a1,a2,a3,b1,b2,x)
df$y = 2*a1+b1+rnorm(n)
```

To perform variable selection, we run the analysis as follows. It will print out the configuration of dap model while running, and the return is a `dap` object.

``` r
test.dap = dap(y~., df)
```

    ## 
    ## ============ DAP Configuration ============
    ## 
    ## INPUT
    ## 
    ##  * individual-level data
    ##  * number of candidate predictors: 1000
    ##  * sample size: 100
    ## 
    ## PROGRAM OPTIONS
    ## 
    ##  * maximum model size allowed [ msize = 1000 ] (no restriction)
    ##  * LD control threshold [ ld_control = 0.25 ]
    ##  * normalizing constant convergence threshold [ converg_thresh = 1.00e-02 ] (log10 scale)
    ##  * number of parallel threads [ thread = 1 ] (OpenMP not available)
    ## 
    ## ===========================================

To see the analysis result by typing the `dap` object name:

``` r
test.dap
```

    ## 
    ## Call:
    ## dap(formula = y ~ ., data = df)
    ## 
    ## Posterior expected model size: 2.025 ( sd = 0.1615 )
    ## LogNC = 63.05 ( Log10NC = 27.38 )
    ## Minimum PIP is estimated at 0.0001001 ( N = 100 )
    ## 
    ## Independent Association Signal Clusters:
    ##    cluster.size  cluster.pip  average.r2  member.predictors
    ## 1             3       0.9999      0.9925           a1 a2 a3
    ## 2             2       0.9999      0.9802              b1 b2
    ## 
    ## One of the best models is:
    ##   y  ~ a1 + b1 
    ## 
    ## Please refer to <dap.object>$signal for PIP of top predictors,
    ##        and <dap.object>$model for configuration of top models.

According to the result, the analysis reveals 2 signal cluster `{a1, a2, a3}` and `{b1, b2}` with high posterior inclusion probability (PIP), as desired. It also suggests one model configuration `y  ~ a1 + b1`. To see information on other plausible model configurations, and details about the predictors within each signal clusters, we can use the `summary` function:

``` r
summary(test.dap)
```

    ## 
    ## Call:
    ## dap(formula = y ~ ., data = df)
    ## 
    ## Posterior expected model size: 2.0251 ( sd = 0.16153 )
    ## LogNC = 63.047 ( Log10NC = 27.381 )
    ## Minimum PIP is estimated at 0.0001001 ( N = 100 )
    ## 
    ## Top Models:
    ##             model  size  posterior   score
    ## 1         a1 + b1     2   0.558895  27.128
    ## 2         a1 + b2     2   0.309274  26.871
    ## 3         a2 + b2     2   0.076997  26.267
    ## 4  a1 + b1 + X969     3   0.015755  25.578
    ## 5         a2 + b1     2   0.015651  25.576
    ## 
    ## 
    ## Independent Association Signal Clusters:
    ##    cluster.size  cluster.pip  r2 matrix          
    ## 1             3      0.99989   0.992522  0.024257
    ## 2             2      0.99989   0.024257  0.980151
    ## 
    ## Cluster 1 :
    ##    predictor       pip  log10abf
    ## 1         a1  0.891022    22.969
    ## 2         a2  0.094350    22.055
    ## 3         a3  0.014513    21.231
    ## 
    ## Cluster 2 :
    ##    predictor      pip  log10abf
    ## 1         b1  0.59091    4.8043
    ## 2         b2  0.40897    4.9174
    ## 
    ## Please refer to <dap.object>$signal for PIP of top predictors,
    ##        and <dap.object>$model for configuration of top models.

### Support of sbams-format file

The sbams format is designed for to represent data from a meta-analytic setting including data from multiple subgroups. The file contains a phenotype section and a genotype section. This package implements an interface `read.sbams` which reads the sbams-format file into a standard R data.frame.

``` r
sbams.file = system.file("sbamsdat", "sim.1.sbams.dat", package = "dap")
sbams.dat  = read.sbams(sbams.file)
```

The data frame generated from `read.sbams` can be seamlessly passed into the `dap` function:

``` r
test.dap.sbams.1 = dap(gene~., sbams.dat, quiet=TRUE)
```

Also, we implement a function `dap.sbams` which can directly read and analyze the sbams-format file:

``` r
test.dap.sbams.2 = dap.sbams(sbams.file, quiet=TRUE)
```

The analysis result `test.dap.sbams.2` is exactly the same as `test.dap.sbams.1`.
