#' TWAS Builder
#'
#' \code{twas} calculates weights on variants.
#'
#' @usage
#' twas(data, dap)
#'
#' @aliases twas
#' @param data  a data frame containing the variables in the model.
#' @param dap   a dap object
#' @return \code{twas} returns an object of \code{"twas"}, which is a list containing the following components:
#' \item{weight}{a data frame of variants and corresponding weights.}
#' \item{ER2}{the expected R-squared.}
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
#' test.dap = dap(y~., df)
#' twas(df, test.dap)
#'
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
twas = function(data, dap){

  # sanity check
  cl = match.call()
  if(class(data)!="data.frame") stop(paste(deparse(cl$data),"is not a data.frame!"))
  if(class(dap) !="dap") stop(paste(deparse(cl$dap)), "is not a dap object!")
  if(dap$model.summary$N != nrow(data))  warning(paste("The number of samples in", deparse(cl$data), "does not match that in", deparse(cl$dap),"!"))
  if(!(dap$model.summary$response %in% names(data))) stop(paste("The response variable", dap$model.summary$response, "is not found in", deparse(cl$data)))

  models = dap$model.summary$model[dap$model.summary$model$size>0,]
  vars = unique(unlist(sapply(models$configuration, function(x) strsplit(x,"\\+"))))
  vars_not_in = !(vars%in% names(data))
  if(sum(vars_not_in) > 0) stop(paste(paste(vars[vars_not_in], collapse=","), "NOT found in", deparse(cl$data)))

  X = data[,vars]
  y = data[,dap$model.summary$response]

  # impute with mean
  y[is.na(y)] = mean(y, na.rm = TRUE)
  X = apply(X, 2, function(t) replace(t, is.na(t), mean(t, na.rm=TRUE)))

  ymean = mean(y)
  const = sum((y-ymean)^2)
  coefs = rep(0, length(vars))
  names(coefs) = vars
  ER2 = 0

  for(i in 1:nrow(models)){
    this_vars = strsplit(models$configuration[i],"\\+")[[1]]
    this_result = lm.fit(cbind(1,as.matrix(X[,this_vars])),y)
    this_weight = models$posterior[i]
    coefs[names(this_result$coef[-1])] = coefs[names(this_result$coef[-1])]+this_result$coef[-1]*this_weight
    ER2 = ER2 + this_weight*(1-sum(this_result$residuals^2)/const)
  }

  index_sorted = sort(abs(coefs), decreasing = TRUE, index.return=TRUE)$ix
  twas_df = data.frame(predictor=names(coefs), weight = coefs)[index_sorted,]
  row.names(twas_df) = 1:length(coefs)
  result = list(coef=twas_df, ER2=ER2)
  class(result) = "twas"
  return(result)
}
