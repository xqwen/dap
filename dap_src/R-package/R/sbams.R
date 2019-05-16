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
#'
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
read.sbams <- function(file) {
  result = .Call(`_dap_read_sbams`, PACKAGE = 'dap', file)
  return(result)
}

