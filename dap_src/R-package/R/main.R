#' Structured Bayesian Model Selection via Deterministic Approximation of Posteriors
#'
#' \code{dap} is used perform Bayesian variable selection among a large scale of predictors, multicolinearty and missingness are allowed. It will automatically impute missingness in the input data with mean values, normalize the response and predictors, and propose top predictive signals and signal clusters according to posterior probability.
#'
#' @usage
#' dap(formula, data, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1, quiet=FALSE)
#'
#' @aliases dap
#' @param formula an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data  a data frame containing the variables in the model.
#' @param ens   (optional) prior expected number of signal clusters, \code{ens=1} by default.
#' @param pi1   (optional) the exchangeable prior probability, values \code{0<pi1<1} accepted. By default -1, \code{pi1=ens/p}, where \code{p} is number of predictors in the input file.
#' @param ld_control (optional) the LD threshold to be considered within a single signal cluster. By default, the threshold is set to 0.25.
#' @param msize (optional) the maximum size of model dap explores. Valid maximum model size ranges from \code{1} to \code{p}. By default -1, it is set to \code{p}, i.e., there is no restriction on how large the true association model can be. If it is specified, the DAP runs DAP-K algorithm and stops at the specified maximum model size.
#' @param converg_thresh (optional)  the stopping condition for model exploration. By default, \code{converg_thresh=0.01}.
#' @param all   (optional) If TRUE, dap will output information for all predictors and all signal clusters. By default, only predictors with \code{PIP > 0.001} and signal clusters with \code{SPIP > 0.25} are output.
#' @param size_limit (optional) the maximum number of predictors allowed in a signal cluster. By default -1, there is no constraint and the size of each signal cluster is completely data determined. Setting a small number forces DAP to cap the number of predictors into each cluster and reduces computation.
#' @param thread (optional) the number of parallel threads to run DAP algorithm, 1 by default. OpenMP is required for multi-thread option.
#' @param quiet (optional) If TRUE, dap will mute running logs.
#' @return \code{dap} returns an object of \code{"dap"}, which is a list containing the following components: \item{cluster}{a data frame, with each line representing the information of one signal cluster, including the size (i.e. number of member predictors), the posterior inclusion probability, and the average LD measures (\eqn{r^2}) for predictors within the cluster.}
#' \item{cluster.r2}{a matrix representing the average LD measures (\eqn{r^2}) for predictors within a cluster and between clusters.}
#' \item{signal}{a data frame of predictors ordered by the posterior inclusion probability (PIP), including predictor name, PIP, log10abf, and the signal cluster it belongs to.}
#' \item{model}{a data frame of models. Specifically, the first column shows the posterior probability of the corresponding model; the second column indicates the size (i.e., the number of predictors) of the model; the third column shows the unnormalized posterior score of the model (defined as \eqn{log10(model prior)+log10(BF)}); and the last column gives the exact configuration of the model.}
#' \item{info}{a list of extra information on the expected model size, sample size and the minimum PIP.}
#' \item{call}{the matched call.}
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
#' test.dap
#'
#' summary(test.dap)
#'
#' @details Please refer to \url{https://github.com/xqwen/dap/tree/master/dap_src} for more details.
#' @seealso \code{\link{summary.dap}} for summaries; and \code{\link{dap.sbams}} for a different interface directly analyzing a sbams-format file.
#' @references Wen, X., Lee, Y., Luca, F., Pique-Regi, R. Efficient Integrative Multi-SNP Association Analysis using Deterministic Approximation of Posteriors. \emph{The American Journal of Human Genetics}, 98(6), 1114--1129
#' @references Lee, Y, Luca, F, Pique-Regi, R,Wen, X. Bayesian Multi-SNP Genetic Association Analysis: Control of FDR and Use of Summary Statistics bioRxiv:316471
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
dap = function(formula, data, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1, quiet=FALSE){
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

  params = list(t=1)
  if(class(ens)=="numeric" & ens > 0)       params$ens=ens
  if(class(pi1)=="numeric" & pi1>0 & pi1<1) params$pi1=pi1
  if(class(ld_control)=="numeric" & ld_control>=0 & ld_control<1) params$ld_control=ld_control
  if(class(msize)=="numeric" & as.integer(msize) >= 1)  params$msize=as.integer(msize)
  if(class(converg_thresh)=="numeric" & converg_thresh>=0) params$converg_thresh=converg_thresh
  if(class(all)=="logical" & all) params$all=1
  if(class(size_limit)=="numeric" & size_limit>=1) params$size_limit=size_limit
  if(class(thread)=="numeric" & as.integer(thread)>1) params$t=as.integer(thread)

  result = .Call(`_dap_dap_sbams`, PACKAGE = 'dap', x, y, 1, params, as.numeric(quiet), all.vars(cl)[1])

  result$call = cl

  class(result) = "dap"
  return(result)
}

#' @export
print.dap = function(object, digits = max(3L, getOption("digits") - 3L)){
  cat("\nCall:\n",
      paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Posterior expected model size:", format(object$info$model.size[1], digits=digits), "( sd =", format(object$info$model.size[2], digits=digits), ")\n")
  cat("LogNC =", format(object$info$log10NC/log10(exp(1)), digits=digits), "( Log10NC =", format(object$info$log10NC, digits=digits), ")\n")
  cat("Minimum PIP is estimated at", format(object$info$PIP.min, digits=digits), "( N =", format(object$info$N, digits=digits), ")\n")

  if(length(object$cluster)){
    cat("\nIndependent Association Signal Clusters:\n")
    signals = sapply(1:nrow(object$cluster), function(t) paste(sort(object$signal[object$signal$cluster==t,"predictor"]), collapse =" "))
    print(format(data.frame(object$cluster, member.predictors=signals), digits=digits), print.gap=2L, quote=FALSE)
  }else{
    cat("\nNo Independent Association Signal Clusters.\n")
  }

  cat("\nOne of the best models is:\n")
  cat("\t", object$info$response,"~", gsub("&", " + ", object$model$configuration[1]), "\n\n")

  cat("Please refer to <dap.object>$signal for PIP of top predictors,\n")
  cat("       and <dap.object>$model for configuration of top models.\n")

  cat("\n")
  invisible(object)
}

#' Summarizing DAP results
#'
#' \code{summary} method for class \code{"dap"}.
#' @aliases summary.dap
#' @aliases print.summary.dap
#' @usage
#' \method{summary}{dap}(object)
#'
#' \method{print}{summary.dap}(object, digits = max(5L, getOption("digits") - 3L))
#' @param object an object of class \code{"dap"}, usually, a result of a call to \code{\link{dap}}.
#' @param digits the number of significant digits to use when printing.
#' @return The function \code{summary.dap} returns a summary of dap results given in \code{object}, using its components \code{"call"} and \code{"info"}, plus \item{top.models}{a data frame for the top 5 (or less if not available) models with highest posterior probabilities, including model configuration, number of predictors involved, the posterior probability of the corresponding model, and the unnormalized posterior score of the model, defined as \eqn{log10(model prior)+log10(BF)}.}
#' \item{clusters}{a data frame for the signal cluster summary, including the size of the cluster (i.e. the number of member predictors), the corresponding posterior inclusion probability, and the remaining columns represent the average LD measures (\eqn{r^2}) for predictors within a cluster and between clusters.}
#' \item{signals}{a list of data frames for predictors in each cluster, including the predictor name, and corresponding posterior inclusion probability.}
#' @export
summary.dap = function(object){
  if(!inherits(object, "dap"))
    warning("calling summary.dap(<fake-dap-object>) ...")
  ans = object[c("call", "info")]
  ans$call = object$call

  n_top_model = min(5, nrow(object$model))
  model = sapply(1:n_top_model, function(i) gsub("&", " + ", object$model$configuration[i]))
  ans$top.models = data.frame(model, object$model[1:n_top_model,c(2,1,3)])

  if(length(object$cluster)){
    ncluster = nrow(object$cluster)
    ans$clusters = data.frame(object$cluster[,-3], object$cluster.r2)
    names(ans$clusters) = c(names(ans$clusters)[1:2], "r2 matrix", rep("", ncluster-1))

    for(i in 1:ncluster){
      ans$signals[[i]] = object$signal[object$signal$cluster==i,-4]
      row.names(ans$signals[[i]]) = 1:nrow(ans$signals[[i]])
    }
  }

  class(ans) = "summary.dap"
  return(ans)
}

#' @export
print.summary.dap = function(object, digits = max(5L, getOption("digits") - 3L)){
  cat("\nCall:\n",
     paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Posterior expected model size:", format(object$info$model.size[1], digits=digits), "( sd =", format(object$info$model.size[2], digits=digits), ")\n")
  cat("LogNC =", format(object$info$log10NC/log10(exp(1)), digits=digits), "( Log10NC =", format(object$info$log10NC, digits=digits), ")\n")
  cat("Minimum PIP is estimated at", format(object$info$PIP.min, digits=digits), "( N =", format(object$info$N, digits=digits), ")\n")

  cat("\nTop Models:\n")
  print(format(object$top.models, digits=digits), print.gap=2L, quote=FALSE)
  cat("\n")


  if(length(object$clusters)){
    cat("\nIndependent Association Signal Clusters:\n")
    print(format(object$clusters, digits=digits), print.gap=2L, quote=FALSE)

    cat("\n")
    for(i in 1:length(object$signals)){
      cat("Cluster",i,":\n")
      print(format(object$signals[[i]], digits=digits), print.gap=2L, quote=FALSE)
      cat("\n")
    }

  }else{
    cat("\nNo Independent Association Signal Clusters.\n")
  }

  cat("Please refer to <dap.object>$signal for PIP of top predictors,\n")
  cat("       and <dap.object>$model for configuration of top models.\n")

  cat("\n")
  invisible(object)
}


#' Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors for SBAMS-format Data
#'
#' \code{dap.sbams} is an interface built especially for SBAMS-format files, which is esigned to perform rigorous enrichment analysis, QTL discovery and multi-SNP fine-mapping analysis in a highly efficient way.
#'
#' @usage
#' dap.sbams(file, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1, quiet=FALSE)
#'
#' @aliases dap.sbams
#' @param file  file path to a sbams file
#' @param ens   (optional) prior expected number of signals, \code{ens=1} by default.
#' @param pi1   (optional) the exchangeable prior probability, values \code{0<pi1<1} accepted. By default -1, \code{pi1=ens/p}, where \code{p} is number of SNPs in the input file.
#' @param ld_control (optional) the LD threshold to be considered within a single signal cluster. By default, the threshold is set to 0.25.
#' @param msize (optional) the maximum size of model dap-g explores. Valid maximum model size ranges from \code{1} to \code{p}. By default -1, it is set to \code{p}, i.e., there is no restriction on how large the true association model can be. If it is specified, the DAP-G runs DAP-K algorithm and stops at the specified maximum model size.
#' @param converg_thresh (optional)  the stopping condition for model exploration. By default, \code{converg_thresh=0.01}.
#' @param all   (optional) If TRUE, dap will output information for all SNPs and all signal clusters. By default, only SNPs with \code{PIP > 0.001} and signal clusters with \code{SPIP > 0.25} are output.
#' @param size_limit (optional) the maximum number of predictors allowed in a signal cluster. By default -1, there is no constraint and the size of each signal cluster is completely data determined. Setting a small number forces DAP to cap the number of predictors into each cluster and reduces computation.
#' @param thread (optional) the number of parallel threads to run DAP algorithm, 1 by default. OpenMP is required for multi-thread option.
#' @param quiet (optional) If TRUE, dap will mute running logs.
#' @return \code{dap} returns an object of \code{"dap"}, which is a list containing the following components: \item{cluster}{a data frame, with each line representing the information of one signal cluster, including the size (i.e. number of member predictors), the posterior inclusion probability, and the average LD measures (\eqn{r^2}) for predictors within the cluster.}
#' \item{cluster.r2}{a matrix representing the average LD measures (\eqn{r^2}) for predictors within a cluster and between clusters.}
#' \item{signal}{a data frame of predictors ordered by the posterior inclusion probability (PIP), including predictor name, PIP, log10abf, and the signal cluster it belongs to.}
#' \item{model}{a data frame of models. Specifically, the first column shows the posterior probability of the corresponding model; the second column indicates the size (i.e., the number of predictors) of the model; the third column shows the unnormalized posterior score of the model (defined as \eqn{log10(model prior)+log10(BF)}); and the last column gives the exact configuration of the model.}
#' \item{info}{a list of extra information on the expected model size, sample size and the minimum PIP.}
#' \item{call}{the matched call.}
#' @details Please refer to \url{https://github.com/xqwen/dap/tree/master/dap_src} for more details.
#' @seealso \code{\link{summary.dap}} for summaries; and \code{\link{read.sbams}} for reading in sbams-format files as an \R \code{data frame} which can call the general version of \code{\link{dap}} function.
#' @references Wen, X., Lee, Y., Luca, F., Pique-Regi, R. Efficient Integrative Multi-SNP Association Analysis using Deterministic Approximation of Posteriors. \emph{The American Journal of Human Genetics}, 98(6), 1114--1129
#' @references Lee, Y, Luca, F, Pique-Regi, R,Wen, X. Bayesian Multi-SNP Genetic Association Analysis: Control of FDR and Use of Summary Statistics bioRxiv:316471
#' @examples \dontrun{
#'
#' sbams.file = system.file("sbamsdat", "sim.1.sbams.dat", package = "dap")
#' test.dap.sbams = dap.sbams(sbams.file)
#' test.dap.sbams
#'
#' summary(test.dap.sbams)
#' }
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
dap.sbams <- function(file, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1, quiet=FALSE)
{
  cl = match.call()
  params = list(data=file)
  if(class(ens)=="numeric" & ens > 0)       params$ens=ens
  if(class(pi1)=="numeric" & pi1>0 & pi1<1) params$pi1=pi1
  if(class(ld_control)=="numeric" & ld_control>=0 & ld_control<1) params$ld_control=ld_control
  if(class(msize)=="numeric" & as.integer(msize) >= 1)  params$msize=as.integer(msize)
  if(class(converg_thresh)=="numeric" & converg_thresh>=0) params$converg_thresh=converg_thresh
  if(class(all)=="logical" & all) params$all=1
  if(class(size_limit)=="numeric" & size_limit>=1) params$size_limit=size_limit
  if(class(thread)=="numeric" & as.integer(thread)>1) params$t=as.integer(thread)

  result = .Call(`_dap_dap_main`, PACKAGE = 'dap', params, as.numeric(quiet))

  result$call = cl
  class(result) = "dap"
  return(result)
}
