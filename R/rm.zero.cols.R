rm.zero.cols = function(xx, eps=10^(-12)) {
  n = ncol(xx)
  vv = NULL
  for(i in 1:n) {
    if(sum(abs(xx[,i])) < eps) {
      vv = c(vv, i)
    }
  }
  if(is.null(vv)) {
    return(xx)
  } else {
    return(xx[,-vv])
  }
}