qvalue.pi0 <-
function(p=NULL, lambda=seq(0.05,0.95,0.05), pi0.method="smoother", smooth.df = 3, smooth.log.pi0 = FALSE, monotone=FALSE, ...) {
	
  if(pi0.method=="smoother") {require(stats)}
  
  ## check arguments
  if(min(p)<0 || max(p)>1) {
    stop("ERROR: p-values not in valid range.")
  }
  if(length(lambda)>1 && length(lambda)<4) {
    stop("ERROR: If length of lambda greater than 1, you need at least 4 values.")
  }
  if(min(lambda) < 0 || max(lambda) >= 1) {
    stop("ERROR: Lambda must be within [0, 1).")
  }
	
  if(length(lambda)==1) {
    pi0 <- mean(p >= lambda)/(1-lambda)
    pi0.lambda = pi0
    pi0 <- min(pi0,1)
  }
  else {
    pi0 <- rep(0,length(lambda))
    for(i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1-lambda[i])
      pi0.lambda = pi0
    }
		if(monotone){ for(x in 2:length(pi0)){ if(pi0[x]>=pi0[x-1] ){ pi0[x] <- pi0[x-1] }} }
		
    if(pi0.method=="smoother") {
      if(smooth.log.pi0) {pi0 <- log(pi0)}
      spi0 <- smooth.spline(lambda,pi0,df=smooth.df)
      pi0 <- predict(spi0,x=max(lambda))$y
      
      if(smooth.log.pi0) {pi0 <- exp(pi0)}
      pi0 <- min(pi0,1)
    }
    else if(pi0.method=="bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0,length(lambda))
      pi0.boot <- rep(0,length(lambda))
      for(i in 1:100) {
        p.boot <- sample(p,size=m,replace=TRUE)
        for(i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot>lambda[i])/(1-lambda[i])
        }
        mse <- mse + (pi0.boot-minpi0)^2
      }
      pi0 <- min(pi0[mse==min(mse)])
      pi0 <- min(pi0,1)
    }
    else {
      stop('ERROR: pi0.method must be one of "smoother" or "bootstrap".')
    }
  }
  if(pi0 <= 0) {
    stop("ERROR: The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda.")
	}
  return(list(pi0=pi0, pi0.lambda=pi0.lambda, lambda=lambda))
}

