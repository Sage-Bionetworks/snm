edge.qvalue <-
function(p,lambda = seq(0, 0.9, 0.05), pi0.method = "smoother",
    fdr.level = NULL, robust = FALSE,smooth.df = 3, smooth.log.pi0 = FALSE, ...) {
  #Function by John Storey for the edge package
	
  err.func <- "edge.qvalue"
  if (min(p,na.rm=TRUE) < 0 || max(p,na.rm=TRUE) > 1) {
    err.msg(err.func,"P-values not in valid range.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && length(lambda) < 4) {
    err.msg(err.func,"If length of lambda greater than 1, you need at least 4 values.")
    return(invisible(1))
  }
  if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 1)) {
    err.msg(err.func,"Lambda must be in [0,1).")
    return(invisible(1))
  }
  m <- length(p)
  if (length(lambda) == 1) {
    if (lambda < 0 || lambda >= 1) {
      err.msg(err.func,"Lambda must be in [0,1).")
      return(invisible(1))
    }
    pi0 <- mean(p >= lambda)/(1 - lambda)
    pi0 <- min(pi0, 1,na.rm=TRUE)
  } else {
    pi0 <- rep(0, length(lambda))
    for (i in 1:length(lambda)) {
      pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
    }
    if (pi0.method == "smoother") {
      if (smooth.log.pi0){ 
        pi0 <- log(pi0)
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
      }
      if (smooth.log.pi0) {
        pi0 <- exp(pi0)
      }
      pi0 <- min(pi0, 1)
    }
    else if (pi0.method == "bootstrap") {
      minpi0 <- min(pi0)
      mse <- rep(0, length(lambda))
      pi0.boot <- rep(0, length(lambda))
      for (i in 1:100) {
        p.boot <- sample(p, size = m, replace = TRUE)
        for (i in 1:length(lambda)) {
          pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - lambda[i])
        }
        mse <- mse + (pi0.boot - minpi0)^2
      }
      pi0 <- min(pi0[mse == min(mse)])
      pi0 <- min(pi0, 1)
    }
    else {
      err.msg(err.func,"'pi0.method' must be one of 'smoother' or 'bootstrap'")
      return(invisible(1))
    }
  }
  if (pi0 <= 0) {
		cat("The estimated pi0 <= 0. Check that you have valid\np-values or use another lambda method.\n")
		#err.msg(err.func,"")
		var<-"error"
		class(var) <- "try-error"
		cat(var,"\n")	
		return(var)
	}
  if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 1)) {
		cat("fdr.level' must be within (0,1]\n")
#		err.msg(err.func,"'fdr.level' must be within (0,1].")
		var<-"error"
		class(var) <- "try-error"
		return(var)
  }
  u <- order(p)
	
  qvalue.rank <- function(x) {
    idx <- sort.list(x)
    fc <- factor(x)
    nl <- length(levels(fc))
    bin <- as.integer(fc)
    tbl <- tabulate(bin)
    cs <- cumsum(tbl)
    tbl <- rep(cs, tbl)
    tbl[idx] <- tbl
    return(tbl)
  }
  v <- qvalue.rank(p)
  qvalue <- pi0 * m * p/v
  if (robust) {
    qvalue <- pi0 * m * p/(v * (1 - (1 - p)^m))
  }
  qvalue[u[m]] <- min(qvalue[u[m]], 1)
  for (i in (m - 1):1) {
    qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]], 1)
  }
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, fdr.level = fdr.level, significant = (qvalue <= fdr.level), lambda = lambda)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, pvalues = p, lambda = lambda)
  }
  class(retval) <- "qvalue"
  return(retval)
}

