sva <-
function(dat, bio.var, adj.var=NULL, n.sv=NULL, num.iter=NULL, diagnose=TRUE, verbose=TRUE) {
  if(class(dat) == "matrix") {
  # It's a matrix which means the user wants to start a new sva session.
    sva.obj = make.sva.obj(dat, bio.var, adj.var, n.sv, num.iter, diagnose)		
	  sva.obj = get.residual.rank(sva.obj)
#		start.id <- sva.obj$iteration
  } else if(class(dat) == "edge") {
  # User didn't define how many iterations they want to update sva fit 
    if(is.null(num.iter)) {
      stop("You are trying to update an existing sva fit.  Please specify the number of iterations.")
    }
    # It's an existing edge object.  Adjust num.iter slot
    dat <- sva.obj$dat
    sva.obj$num.iter <- sva.obj$iteration + sva.obj$num.iter
    singular.values = matrix(NA, nr=sva.obj$num.iter, nc=ncol(sva.obj$dat))		
    singular.values[1:nrow(sva.obj$singular.values),] <- sva.obj$singular.values 
    sva.obj$singular.values <- singular.values
  } else if(class(dat) != "edge") { 
    # User variable is not valid
    stop("dat must be a matrix or an object of class sva.")
  }
	start.id <- (sva.obj$iteration+1)
  for(i in start.id:sva.obj$num.iter) {
    if(verbose){cat("\r", "Iteration: ", i)}		
    # Estimate Pr(b_i = 0)		
    mod.b <- sva.obj$bio.var
    mod0.b <- cbind(sva.obj$adj.var,sva.obj$svd[[i-1]]$v)
		sva.obj.b = make.sva.obj(dat, mod.b, mod0.b, n.sv, num.iter, diagnose)		
    fit <- edge.fit(sva.obj.b)
    ptmp.b <- edge.glr(fit,sva.obj.b$df.full,sva.obj.b$df.null,norm.pval=TRUE)$pval
    pprob.b <- (1-qvalue.lfdr(ptmp.b))
		
    # Estimate Pr(u_i = 0)	
    mod.u <- model.matrix(~sva.obj$svd[[i-1]]$v)
    mod0.u <- sva.obj$adj.var
    sva.obj.u = make.sva.obj(dat, mod.u, mod0.u, n.sv, num.iter, diagnose)		
    fit <- edge.fit(sva.obj.u)
    ptmp.u <- edge.glr(fit,sva.obj.u$df.full,sva.obj.u$df.null,norm.pval=TRUE)$pval
    pprob.u <- (1-qvalue.lfdr(ptmp.u))

    # Estimate Pr(u_i = 0) * (1 - Pr(b_i = 0))
    pprob <- pprob.u * (1 - pprob.b);
		
    # Perform weighted SVD.  If no adjustment variables decompose weighted data
    if(ncol(sva.obj$adj.var) == 1) { 
      B = crossprod(sweep(pprob*sva.obj$dat, 1, rowMeans(pprob*sva.obj$dat)))
    }else{
      # If at least one adjustment variable decompose weighted residuals obtained
      # by removing adjustment variables.
      mod <- sva.obj$adj.var
      H <- mod %*% solve(t(mod) %*% mod) %*% t(mod)
      res <- sva.obj$dat - t(H %*% t(sva.obj$dat))			
      B = crossprod(pprob*res)			
    }
    s = svd(B, nu = 0,nv=sva.obj$n.sv)

    # Update sva.obj 
    ds <- s$d / sum(s$d)
    sva.obj$singular.values[i,] <- ds
    sva.obj$svd[[i]] <- s
    sva.obj$iteration <- sva.obj$iteration + 1
    if(sva.obj$diagnose) { sva.diagnostic.plot(ptmp.b, pprob.b, ptmp.u, pprob.u, pprob, sva.obj)}
  }
	cat("\n")
  sva.obj
}

