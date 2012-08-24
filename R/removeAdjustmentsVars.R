remove.adjustment.vars <- function(snm.obj) {
  dat = snm.obj$dat
  n = snm.obj$n.arrays
  xx0 = snm.obj$adj.var
  df0 = snm.obj$df.null
  xx1 = cbind(snm.obj$adj.var, snm.obj$bio.var)
  P1 = solve(t(xx1) %*% xx1) %*% t(xx1)
  cfs1 = t(P1 %*% t(dat))
  res1 = dat - cfs1 %*% t(xx1)
  x1 <- cbind(rep(1,dim(dat)[2]), snm.obj$bio.var)
  if(ncol(snm.obj$adj.var) > 1) {
    cfs2 <- cfs1[, -(2:ncol(snm.obj$adj.var))]
  }else{
    cfs2 <- cfs1
  }
  fit1 <- cfs2 %*% t(x1)
  snm.obj$dat <- fit1 + res1
  snm.obj
}

