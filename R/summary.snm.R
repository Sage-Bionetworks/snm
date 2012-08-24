snm.summary <- function(object, cuts=c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1), ...) {

  cat('\n'); cat('SNM Data and Model Summary', '\n', '\n')
  cat('Total number of arrays:', ncol(object$norm.dat), '\n')
  cat('Total number of probes:', nrow(object$norm.dat), '\n', '\n')
  if(!is.null(object$bio.var)) {
    cat('Final estimated proportion of null probes: '); cat(round(object$pi0, 3)); cat('\n'); cat('\n') 
    cat("Cumulative number of significant calls:\n")
    qobj = edge.qvalue(object$pval, ...)
    counts = sapply(cuts, function(x) c("p-value"=sum(qobj$pvalues < x), "q-value"=sum(qobj$qvalues < x)))
    colnames(counts) = paste("<", cuts, sep="")
    print(counts)
    cat("\n")
    cat('Full model degrees of freedom:', object$df1, '\n')
    cat('Null model degrees of freedom:', object$df0, '\n', '\n')
    cat('Biological variables:', '\n'); print(signif(t(object$bio.var)), digits=3); cat('\n')
  }
  cat('Adjustment variables:', '\n'); print(signif(t(object$adj.var)), digits=3); cat('\n') 
  cat('Intensity-dependent variables:', '\n'); print(t(object$int.var)); cat('\n')
  cat('Function call:', '\n'); print(object$call)
  
}

summary.snm <- function(object, ...) {
  snm.summary(object, ...)
}
