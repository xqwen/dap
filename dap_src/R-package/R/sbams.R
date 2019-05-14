#' Read a sbams file
#'
#' This function will import the sbams file, automatically impute missing values with mean, and regress pheno and geno out of controlled covariants if applicable.
#'
#' @param file file path to the sbams file
#' @return a sbams object with pheno and geno, and controlled covariates if applicable.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
read.sbams <- function(file) {
  result = .Call(`_dap_read_sbams`, PACKAGE = 'dap', file)
  class(result) = "sbams"
  return(result)
}
