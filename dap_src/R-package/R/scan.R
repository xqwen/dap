#' Single-SNP Association Analysis
#'
#' \code{scan} is used perform perform single-SNP fine-mapping analysis.
#'
#' @usage
#' scan(formula, data)
#'
#' @aliases scan
#' @param formula an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data  a data frame containing the variables in the model.
#' @return \code{scan} returns a \code{data.frame} containing effect estimates for each single SNP.
#' @examples
#' set.seed(0)
#' n = 100
#' p = 1000
#'
#' a1 = rnorm(n)
#' a2 = 10*a1+rnorm(n)
#' a3 = 2*a1+9*a2+rnorm(n)
#'
#' b1 = rnorm(n)
#' b2 = 8*b1+rnorm(n)
#'
#' x = matrix(rnorm(n*(p-5)), nrow=n)
#'
#' df = data.frame(a1,a2,a3,b1,b2,x)
#' df$y = 2*a1+b1+rnorm(n)
#'
#' test.scan = scan(y~., df)
#' head(test.scan)
#'
#' @seealso \code{\link{scan.sbams}} for a different interface directly analyzing a sbams-format file.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
scan = function(formula, data){
  cl = match.call()
  mf = match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")

  y = model.response(mf, "numeric")
  attr(mt, "intercept") = 0
  x <- model.matrix(mt, mf)

  # impute with mean
  y[is.na(y)] = mean(y, na.rm = TRUE)
  x = apply(x, 2, function(t) replace(t, is.na(t), mean(t, na.rm=TRUE)))

  params = list(x=x,y=y,pheno_name=all.vars(cl)[1],scan=TRUE)

  result = .Call(`_dap_dap_main`, PACKAGE = 'dap', 2, params, 1)

  return(result[[1]])
}

#' Single-SNP Association Analysis for SBAMS-format Data
#'
#' \code{scan.sbams} is an interface built especially for SBAMS-format files, which is designed to perform single-SNP fine-mapping analysis.
#'
#' @usage
#' scan.sbams(file)
#'
#' @aliases scan.sbams
#' @param file  file path to a sbams file
#' @return \code{scan.sbams} returns a \code{data.frame} containing effect estimates for each single SNP.
#' @examples \dontrun{
#'
#' sbams.file = system.file("sbamsdat", "sim.1.sbams.dat", package = "dap")
#' test.scan.sbams = scan.sbams(sbams.file)
#' head(test.scan.sbams)
#'
#' }
#' @seealso \code{\link{read.sbams}} for reading in sbams-format files as an \R \code{data.frame} which can call the general version of \code{\link{scan}} function.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
scan.sbams <- function(file)
{
  params = list(data=file, scan=TRUE)

  result = .Call(`_dap_dap_main`, PACKAGE = 'dap', 1, params, 1)

  return(result[[1]])
}
