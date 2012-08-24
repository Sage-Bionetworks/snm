snm <-
		function(raw.dat, bio.var=NULL, adj.var=NULL, int.var=NULL,
        weights=NULL, spline.dim = 4, num.iter = 10, nbins=50,
        rm.adj=FALSE, verbose=TRUE, diagnose=TRUE)
{
  if(is.null(bio.var)) {
    warning("bio.var=NULL, so all probes will be treated as 'null' in the normalization.")
  }
  
  snm.obj <- make.snm.obj(Y=raw.dat, bio.var, adj.var, int.var, spline.dim,
      nbins, weights, diagnose, rm.adj)  
  if(!is.null(snm.obj$int.var)) {		  
    #Set up model fitting loop
    basisSplineFunction <- buildBasisFunction(snm.obj)
    pi0s <- rep(0,num.iter)
    obs.fit <- edge.fit(snm.obj, odp=FALSE)
		if(is.null(snm.obj$bio.var)) { 
      snm.obj$pi0 <- 1
      snm.obj$nulls <- 1:snm.obj$n.probes
      obs.fit = list()
      xx0 = snm.obj$adj.var
      H0 = xx0 %*% solve(t(xx0) %*% xx0) %*% t(xx0)
      obs.fit$fit0 = t(H0 %*% t(snm.obj$dat))
      obs.fit$res1 = snm.obj$dat - obs.fit$fit0
      snm.obj <- suppressWarnings(fit.model(obs.fit, snm.obj, basisSplineFunction))
      snm.obj$dat <- snm.obj$r.dat - snm.obj$array.fx
      snm.obj$pi0s <- NULL
    } else {
    	#Otherwise perform model fitting loop
      pi0s <- rep(0,num.iter)
      obs.fit <- edge.fit(snm.obj, odp=FALSE)
      for (i in 1:num.iter) {
        if(verbose){cat("\r", "Iteration: ", i)}
        obs.stat <- edge.glr(obs.fit, df1=snm.obj$df.full, df0=snm.obj$df.null, norm.pval=TRUE)
      	snm.obj$pi0 <- edge.qvalue(obs.stat$pval)$pi0
      	snm.obj$nulls <- calculate.nulls(obs.stat$pval, snm.obj$pi0)
        snm.obj <- suppressWarnings(fit.model(obs.fit, snm.obj, basisSplineFunction))
        snm.obj$dat <- snm.obj$r.dat - snm.obj$array.fx
      	obs.fit <- edge.fit(snm.obj, odp=FALSE)
      	if(diagnose) {
          snm.diagnostic.plot(obs.fit,obs.stat$pval,snm.obj,iter=i)
      	}
      	pi0s[i] <- snm.obj$pi0
      }
      snm.obj$pi0s <- pi0s
    }
    #Extract subset of array effects (to reduce size of return object)
    #Set seed before running snm if this needs to be reproduced
    if(dim(snm.obj$dat)[1] > (nbins*100)) {
      th <- sample(dim(snm.obj$dat)[1], (nbins*100))
    } else {
      th <- 1:dim(snm.obj$dat)[1]
    } 
    M.ret = snm.obj$M[th,]
    array.fx.ret = snm.obj$array.fx[th,]
    for(i in 1:ncol(M.ret)) {
      oo <- order(M.ret[,i])
      M.ret[,i] = M.ret[oo,i]
      array.fx.ret[,i] = array.fx.ret[oo,i]
    }
    #Get final model fit significance
    if(is.null(bio.var)) { 
      obs.stat <- NULL
      snm.obj$pval <- NULL
      snm.obj$pi0 <- NULL
    } else {
      obs.stat <- edge.glr(obs.fit, df1=snm.obj$df.full, df0=snm.obj$df.null, norm.pval=TRUE)
      snm.obj$pval <- obs.stat$pval
      snm.obj$pi0 <- edge.qvalue(obs.stat$pval)$pi0
      if(snm.obj$rm.adj){
      	snm.obj <- remove.adjustment.vars(snm.obj)
      }
    }
  } else {
    #No intensity-dep effects
		if(!is.null(snm.obj$bio.var)) { 
	    obs.fit <- edge.fit(snm.obj, odp=FALSE)
  	  obs.stat <- edge.glr(obs.fit, df1=snm.obj$df.full, df0=snm.obj$df.null, norm.pval=TRUE)
    	if(!snm.obj$rm.adj) {
      	stop("int.var=NULL and rm.adj=FALSE, so there is nothing to do.")
    	}
    	snm.obj <- remove.adjustment.vars(snm.obj)
    	snm.obj$pval <- obs.stat$pval; 
    	snm.obj$pi0s <- edge.qvalue(obs.stat$pval)$pi0
    	snm.obj$pi0 <- snm.obj$pi0s
		}else{
			# No biological effects, so just remove the adjustment variable
			snm.obj <- remove.adjustment.vars(snm.obj)
			snm.obj$pval <- NULL
			snm.obj$pi0s <- NULL
			snm.obj$pi0 <- NULL
		}
   	M.ret = NULL
    array.fx.ret = NULL
  }
	
  dimnames(snm.obj$adj.var) = snm.obj$dimnames.adj
  dimnames(snm.obj$bio.var) = snm.obj$dimnames.bio
  
  snm.ret <- list(norm.dat=snm.obj$dat, pval=snm.obj$pval, pi0=snm.obj$pi0, iter.pi0s=snm.obj$pi0s, nulls=snm.obj$nulls,
      M=M.ret, array.fx=array.fx.ret, bio.var=snm.obj$bio.var, adj.var=snm.obj$adj.var, int.var=snm.obj$int.var,
      df1=snm.obj$df.full, df0=snm.obj$df.null, raw.dat=raw.dat, rm.adj=rm.adj, call = match.call())
	
  class(snm.ret) = "snm"
	if(verbose){cat("\n")}
  return(snm.ret)
}