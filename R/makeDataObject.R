makeDataObject <-
		function(Y.pooled, np, snm.obj, exp,bins) 
{
  lnp <- length(np)
  num.arrays <- dim(Y.pooled)[2]
  signal <- data.frame(y = as.numeric(Y.pooled), 
      probes = factor(rep(1:lnp, times = num.arrays)))
  f <- lapply(1:length(snm.obj$int.var), function(x) {
    		i <- snm.obj$int.var[, x]
    		rep(i, each = lnp)
  		})
  f <- as.data.frame(f)
  colnames(f) <- names(snm.obj$int.var)
  if (ncol(snm.obj$adj.var) > 1) {
    z <- lapply(2:(ncol(snm.obj$adj.var)), function(x) {
      		i <- snm.obj$adj.var[, x]
      		rep(i, each = lnp)
    		})
    z <- as.data.frame(z)
    colnames(z) <- colnames(snm.obj$adj.var)[-1]
    df <- cbind(cbind(signal, z), f)
  }else {
    df <- cbind(signal, f)
  }
  df$weights = sapply(bins,length) / 50 #mean(sapply(bins,length))
  df$weights[df$weights >1] <- 1
  TMP <- sapply(1:dim(df)[2], function(x) {
    		assign(colnames(df)[x], df[, x], envir = exp)
  		})
  as.data.frame(df)
}
