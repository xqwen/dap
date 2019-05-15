#' Read a sbams file
#'
#' This function will import the sbams file, automatically impute missing values with mean, and regress pheno and geno out of controlled covariants if applicable.
#'
#' @param file file path to the sbams file
#' @return a data.frame with the first colunm as the normalized phenotype, and the following columns as the normalized genotypes.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
read.sbams <- function(file) {
  result = .Call(`_dap_read_sbams`, PACKAGE = 'dap', file)
  return(result)
}

