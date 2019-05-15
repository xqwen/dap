#' dap
#'
#' Structured Bayesian Model Selection via Deterministic Approximation of Posteriors
#'
#' @usage
#' dap(formula, data, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1)
#'
#' @param formula an object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data  a data frame containing the variables in the model.
#' @param ens   (optional) prior expected number of signal clusters, \code{ens=1} by default.
#' @param pi1   (optional) the exchangeable prior probability, values \code{0<pi1<1} accepted. By default -1, \code{pi1=ens/p}, where \code{p} is number of predictors in the input file.
#' @param ld_control (optional) the LD threshold to be considered within a single signal cluster. By default, the threshold is set to 0.25.
#' @param msize (optional) the maximum size of model dap explores. Valid maximum model size ranges from \code{1} to \code{p}. By default -1, it is set to \code{p}, i.e., there is no restriction on how large the true association model can be. If it is specified, the DAP runs DAP-K algorithm and stops at the specified maximum model size.
#' @param converg_thresh (optional)  the stopping condition for model exploration. By default, \code{converg_thresh=0.01}.
#' @param all   (optional) If TRUE, dap will output information for all predictors and all signal clusters. By default, only predictors with \code{PIP > 0.001} and signal clusters with \code{SPIP > 0.25} are output.
#' @param thread (optional) the number of parallel threads to run DAP algorithm, 1 by default. OpenMP is required for multi-thread option.
#' @return a dap object including model summary, SNP summary and signal cluster summary
#' @details Please refer to \url{https://github.com/xqwen/dap/tree/master/dap_src} for more details.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
dap = function(formula, data, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1){
  cl = match.call()
  mf = match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  y = model.response(mf, "numeric")
  mt <- attr(mf, "terms")
  attr(mt, "intercept") = 0
  x <- model.matrix(mt, mf, contrasts)

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

  result = .Call(`_dap_dap_sbams`, PACKAGE = 'dap', x, y, 1, params)

  result$call = cl

  class(result) = "dap"
  return(result)
}

print.dap = function(object, digits = max(3L, getOption("digits") - 3L)){
  cat("\nCall:\n",
      paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Posterior expected model size:", format(object$info$model.size[1], digits=digits), "( sd =", format(object$info$model.size[2], digits=digits), ")\n")
  cat("LogNC =", format(object$info$log10NC/log10(exp(1)), digits=digits), "( Log10NC =", format(object$info$log10NC, digits=digits), ")\n")
  cat("Minimum PIP is estimated at", format(object$info$PIP.min, digits=digits), "( N =", format(object$info$N, digits=digits), ")\n")

  if(length(object$cluster)){
    cat("\nIndependent Association Signal Clusters:\n")
    signals = sapply(1:nrow(object$cluster), function(t) paste(sort(object$signal[test1$signal$cluster==t,"predictor"]), collapse =" "))
    print(format(data.frame(object$cluster, member.predictors=signals), digits=digits), print.gap=2L, quote=FALSE)
  }else{
    cat("\nNo Independent Association Signal Clusters.\n")
  }

  cat("\n")
  invisible(object)
}

summary.dap = function(object){
  if(!inherits(object, "dap"))
    warning("calling summary.dap(<fake-dap-object>) ...")
  ans = object[c("call", "info")]
  ans$call = object$call

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

print.summary.dap = function(object, digits = max(5L, getOption("digits") - 3L)){
  cat("\nCall:\n",
     paste(deparse(object$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Posterior expected model size:", format(object$info$model.size[1], digits=digits), "( sd =", format(object$info$model.size[2], digits=digits), ")\n")
  cat("LogNC =", format(object$info$log10NC/log10(exp(1)), digits=digits), "( Log10NC =", format(object$info$log10NC, digits=digits), ")\n")
  cat("Minimum PIP is estimated at", format(object$info$PIP.min, digits=digits), "( N =", format(object$info$N, digits=digits), ")\n")

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

  cat("Please refer to <dap.object>$signal for pip of top predictors,\n")
  cat("       and <dap.object>$model for configuration of top models.\n")

  cat("\n")
  invisible(object)
}


#' dap.sbams
#'
#' Integrative Genetic Association Analysis using Deterministic Approximation of Posteriors for SBAMS-format Data
#'
#' @usage
#' dap.sbams(file, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1)
#'
#' @param file  file path to a sbams file
#' @param ens   (optional) prior expected number of signals, \code{ens=1} by default.
#' @param pi1   (optional) the exchangeable prior probability, values \code{0<pi1<1} accepted. By default -1, \code{pi1=ens/p}, where \code{p} is number of SNPs in the input file.
#' @param ld_control (optional) the LD threshold to be considered within a single signal cluster. By default, the threshold is set to 0.25.
#' @param msize (optional) the maximum size of model dap-g explores. Valid maximum model size ranges from \code{1} to \code{p}. By default -1, it is set to \code{p}, i.e., there is no restriction on how large the true association model can be. If it is specified, the DAP-G runs DAP-K algorithm and stops at the specified maximum model size.
#' @param converg_thresh (optional)  the stopping condition for model exploration. By default, \code{converg_thresh=0.01}.
#' @param all   (optional) If TRUE, dap will output information for all SNPs and all signal clusters. By default, only SNPs with \code{PIP > 0.001} and signal clusters with \code{SPIP > 0.25} are output.
#' @param thread (optional) the number of parallel threads to run DAP algorithm, 1 by default. OpenMP is required for multi-thread option.
#' @return a dap object including model summary, SNP summary and signal cluster summary
#' @details Please refer to \url{https://github.com/xqwen/dap/tree/master/dap_src} for more details.
#' @useDynLib dap, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
dap.sbams <- function(file, ens=1, pi1=-1, ld_control=0.25, msize=-1, converg_thresh=0.01, all=FALSE, size_limit=-1, thread=1)
{
  params = list(data=file)
  if(class(ens)=="numeric" & ens > 0)       params$ens=ens
  if(class(pi1)=="numeric" & pi1>0 & pi1<1) params$pi1=pi1
  if(class(ld_control)=="numeric" & ld_control>=0 & ld_control<1) params$ld_control=ld_control
  if(class(msize)=="numeric" & as.integer(msize) >= 1)  params$msize=as.integer(msize)
  if(class(converg_thresh)=="numeric" & converg_thresh>=0) params$converg_thresh=converg_thresh
  if(class(all)=="logical" & all) params$all=1
  if(class(size_limit)=="numeric" & size_limit>=1) params$size_limit=size_limit
  if(class(thread)=="numeric" & as.integer(thread)>1) params$t=as.integer(thread)

  result = .Call(`_dap_dap_main`, PACKAGE = 'dap', params)
  cl = match.call()
  result$call = cl
  class(result) = "dap"
  return(result)
}
