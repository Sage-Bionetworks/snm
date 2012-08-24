qvalue.lfdr <-
function(p, pi0=NULL, trunc=TRUE, monotone=TRUE, transf=c("probit", "logit"), adj=1.5, eps=10^-8, ...) {
	
  require(stats)
  
  ## check arguments
  if(min(p)<0 || max(p)>1) {
    print("ERROR: p-values not in valid range.")
    return(0)
  }
  
  if(is.null(pi0)) pi0 = qvalue.pi0(p, ...)$pi0
  n = length(p)
  transf = match.arg(transf)
  
  if(transf=="probit") {
    p = pmax(p, eps)
    p = pmin(p, 1-eps)
    x = qnorm(p)
    myd = density(x, adjust=adj)
    mys = smooth.spline(x=myd$x, y=myd$y)
    y = predict(mys, x)$y
    lfdr = pi0*dnorm(x)/y
  }
  
  if(transf=="logit") {
    x = log((p+eps)/(1-p+eps))
    myd = density(x, adjust=adj)
    mys = smooth.spline(x=myd$x, y=myd$y)
    y = predict(mys, x)$y
    dx = exp(x)/(1+exp(x))^2
    lfdr = pi0*dx/y
  }
  
  if(trunc) {lfdr[lfdr > 1] = 1}
  if(monotone) {	
    lfdr = lfdr[order(p)]
    for(i in 2:n) {
      if(lfdr[i] < lfdr[i-1]) {lfdr[i] = lfdr[i-1]}
    }
    lfdr = lfdr[rank(p)]
  }
  
  return(lfdr)
}

