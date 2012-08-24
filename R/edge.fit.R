edge.fit <-
function(edge.obj,odp=FALSE) {
  #Function by John Storey for the edge package
  
  err.func = "edge.fit"
  
  dat = edge.obj$dat
  n = edge.obj$n.arrays
  xx0 = edge.obj$adj.var
  df0 = edge.obj$df.null
  if(!odp) {
    xx1 = cbind(edge.obj$adj.var, edge.obj$bio.var)
  } else {
    xx1 = edge.obj$bio.var
  }
  df1 = edge.obj$df.full
  
  if(!is.null(edge.obj$ind)) {
    xxi = edge.obj$ind
    Hi = xxi %*% solve(t(xxi) %*% xxi) %*% t(xxi)
    fit.ind = t(Hi %*% t(dat))
    dat = dat - fit.ind
    xx1 = xx1 - Hi %*% xx1
    xx0 = xx0 - Hi %*% xx0
    xx1 = rm.zero.cols(xx1)
    xx0 = rm.zero.cols(xx0)
  }
  
  H0 = xx0 %*% solve(t(xx0) %*% xx0) %*% t(xx0)
  fit0 = t(H0 %*% t(dat))
  res0 = dat - fit0
  
  if(odp) {
    xx1 = xx1 - H0 %*% xx1
    H1 = xx1 %*% solve(t(xx1) %*% xx1) %*% t(xx1)
    fit1 = t(H1 %*% t(res0))
    res1 = res0-fit1
    var1 = rowSums(res1^2)/(n-df1)
    var0 = rowSums(res0^2)/(n-df0)
    return(list(fit1=fit1, fit0=fit0, res1=res1, res0=res0, var1=var1, var0=var0, ind=edge.obj$ind, odp=odp))
  } else {
    H1 = xx1 %*% solve(t(xx1) %*% xx1) %*% t(xx1)
    C1 = solve(t(xx1) %*% xx1) %*% t(xx1) %*% t(dat)
    fit1 = t(H1 %*% t(dat))
    res1 = dat - fit1
    return(list(fit0=fit0, res1=res1, res0=res0, ind=edge.obj$ind, odp=odp, c1=C1))
  }
}

