#' Read a sbams file
#'
#' This function will import the sbams file, automatically impute missing values with mean, regress pheno and geno out of controlled covariants if applicable, and generate an \R data frame which can be passed into the \code{\link{dap}} function.
#'
#' @param file file path to the sbams file
#' @return a data.frame with the first colunm as the normalized phenotype, and the following columns as the normalized genotypes.
#' @details Please refer to \url{https://github.com/xqwen/dap/tree/master/dap_src} for more details.
#' @examples
#' sbams.file = system.file("sbamsdat", "sim.1.sbams.dat", package = "dap")
#' sbams.dat  = read.sbams(sbams.file)
#' test.dap.sbams.dat = dap(gene~., sbams.dat)
#'
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
read.sbams <- function(file) {
  result = .Call(`_dap_read_sbams`, PACKAGE = 'dap', file)
  return(result)
}

#' Extract Sufficient Statistics from a sbams file
#' 
#' The sufficient summary statistics refer to the following information: \itemize{
#' \item estimated effect size and corresponding estimated standard error from single SNP testing for each SNP
#' \item correlation matrix of SNPs
#' \item sample size
#' \item total sum of squares for the quantitative phenotype
#' }
#' 
#' @param file file path to the sbams file
#' @param ens   (optional) prior expected number of signals, \code{ens=1} by default.
#' @param pi1   (optional) the exchangeable prior probability, values \code{0<pi1<1} accepted. By default -1, \code{pi1=ens/p}, where \code{p} is number of SNPs in the input file.
#' @return \code{extract.sbams} returns a list containing the following components: \item{est}{a \code{data.frame} including snp names, estimated effect size and corresponding estimated standard error.}
#' \item{ld}{a \code{matrix} representing the correlation matrix of SNPs.}
#' \item{n}{sample size.}
#' \item{syy}{total sum of squares for the quantitative phenotype.}
#' \item{name}{name of the phenotype}
#' @examples
#' sbams.file = system.file("sbamsdat", "sim.1.sbams.dat", package = "dap")
#' ss = extract.sbams(sbams.file)
#' dap.ss(ss$est, ss$ld, ss$n, ss$syy, ss$name)
#' 
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
extract.sbams <- function(file, ens=1, pi1=-1) {
  params = list(data=file)
  
  if(class(ens)=="numeric" & ens > 0)       params$ens=ens
  if(class(pi1)=="numeric" & pi1>0 & pi1<1) params$pi1=pi1
  
  result = .Call(`_dap_extract_sbams`, PACKAGE = 'dap', params)
  return(result)
}

