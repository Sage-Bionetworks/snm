snm.plot <- function(x, col.by=NULL, ...) {

  #Process col.by argument, which is used for array effects only
  if(is.null(col.by)) {
    col.by = rep(1,ncol(x$M))
  } else {
    if(is.matrix(col.by)) {
      if(ncol(col.by)==1) {
        col.by = as.vector(col.by)
      } else if(!all(as.vector(col.by)==0 + (as.vector(col.by)==1))) {
        stop("If col.by is a matrix, then it can only take values 0 or 1.")
      } else {
        col.by = drop(col.by %*% (1:ncol(col.by)))
      }
    }
    col.by = as.factor(col.by)
    k = length(levels(col.by))
    col.by = drop(model.matrix(~-1+col.by) %*% (1:k))
    if(k < 4) {
      col.by = c("red", "blue", "green")[col.by]
    } else {
      col.by = rainbow(k)[col.by]
    }
  }
  
  if(is.null(x$bio.var) && !is.null(x$int.var)) {
    plot(x$M[,1], x$array.fx[,1], type="l", lwd=1, xlim=range(x$M), ylim=range(x$array.fx), col=col.by[1],
         xlab="Estimated Intensity", ylab="Estimated Effect by Array",main="Intensity-Dependent Effects")
    for(i in 2:dim(x$M)[2]) {
      points(x$M[,i], x$array.fx[,i],type="l", col=col.by[i])
    }
  } else {
  
    par(mfrow=c(2,2), oma=c(1,0,2,0))
    
  #Plot A
    pi0.diffs = x$iter.pi0s - x$pi0 
    plot(pi0.diffs, pch=19, ylim=c(-max(abs(pi0.diffs)), max(abs(pi0.diffs))),
         xlab="Iteration", ylab="pi0.iteration - pi0.final", main="Convergence")
    abline(h=0, lty=2)

  #Plot B
    if(x$rm.adj) {
      plot(0, 0, xlab=" ", ylab=" ", main="Latent Structure", pch=" ")
      warning("It is not possible to plot latent structure when rm.adj=TRUE.")
    } else {
      res1 <- x$norm.dat - fitted(x)$fit1
      u <- fast.svd(res1,tol=0)
      plot(100*u$d^2 / sum(u$d^2), ylim=c(0,100), pch=19,
           main="Residual Latent Structure",
           xlab="Principal Component",
           ylab="Percent Variation Explained")
    }
        
  #Plot C
    if(is.null(x$M)) {
      plot(0, xlab=" ", ylab=" ", main="Intensity-Dependent Effects")
      warning("No intensity-dependent effects were estimated.")
    } else {
      plot(x$M[,1], x$array.fx[,1],
           type="l", lwd=1, xlim=range(x$M), ylim=range(x$array.fx), col=col.by[1],
           xlab="Estimated Intensity", ylab="Estimated Effect by Array",main="Intensity-Dependent Effects")
      for(i in 2:dim(x$M)[2]) {
        points(x$M[,i], x$array.fx[,i],type="l", col=col.by[i])
      }
    }
  
  #Plot D
    hist(x$pval, xlab="P-values",main="P-value Distribution", freq=FALSE)
    abline(v=min(x$pval[x$nulls]),col="red")
    pi0 = round(x$pi0,3)
    abline(h=pi0, lty=3)
    mtext(substitute(hat(pi)[0] == that, list(that=pi0)))

    #Overall Title
    title("SNM Diagnostic Plot", outer=TRUE)
  }

}

plot.snm <- function(x, ...) {
  snm.plot(x, ...)
}
