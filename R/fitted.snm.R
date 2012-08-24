snm.fitted <- function(object, ...) {

  if(object$rm.adj) {
    stop("This function is not permitted when rm.adj=TRUE.")
  }

  if(is.null(object$bio.var)) {
    stop("This function is not permitted when bio.var=NULL.")
  }
  
  edge.obj = list(dat=object$norm.dat, bio.var=object$bio.var,
    adj.var=object$adj.var, df.full=object$df1, df.null=object$df0,
    n.arrays=ncol(object$norm.dat), ind=NULL)
  class(edge.obj) = "edge"
  snm.fit = edge.fit(edge.obj, odp=FALSE)
  fit1 = object$norm.dat - snm.fit$res1
  fit0 = snm.fit$fit0

  ret.val = list(fit0=fit0, fit1=fit1)
  return(ret.val)

}

fitted.snm <- function(object, ...) {
  snm.fitted(object, ...)
}

