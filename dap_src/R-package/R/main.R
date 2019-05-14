#' Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors
#'
#' The \code{dap} function accepts two different sets of inputs.
#' \enumerate{
#' \item sbams (either file path or object)
#' \item summary statistics (ld, est, n, syy)
#' }
#' @usage
#' dap(file="/path/to/sbams/file")
#' dap(sbams=a.sbams.object)
#' dap(ld="/path/to/ld/file", est="/path/to/est/file", n, syy)
#'
#' @param file  file path to a sbams file
#' @param sbams a sbams object
#' @param ld    file path to a ld file
#' @param est   file path to a est file
#' @param n     an integer indicating the sample size
#' @param syy   a float number indicating the total sum of squares for the quantitative trait
#' @return a dap object including model summary, SNP summary and signal cluster summary
#'
#'
#' @details The summary statistics can be obtained by \code{summary(a.sbams.object)}.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
dap <- function(file=NULL, sbams=NULL, ld=NULL, est=NULL, n=NA, syy=NA)
{
  # if(!is.null(sbams))
  # {
  #   stopifnot(class(sbams)=="sbams")
  #   return(.Call(`_dap_dap_main`, PACKAGE = 'dap', list(data=sbams$file)))
  # }
  if(!is.null(file))
  {
    result = .Call(`_dap_dap_main`, PACKAGE = 'dap', list(data=file))
    class(result) = "dap"
    return(result)
  }
  # if(!is.null(ld) & !is.null(est) & !is.na(n) & !is.na(syy))
  # {
  #   return(.Call(`_dap_dap_main`, PACKAGE = 'dap', list(ld=ld,est=est,n=n,syy=syy)))
  # }
}
