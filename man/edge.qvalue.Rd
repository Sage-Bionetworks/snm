\name{edge.qvalue}
\alias{edge.qvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Qvalue function from EDGE software
}
\description{
Implement Qvalue methodology on list of P values.
}
\usage{
edge.qvalue(p, lambda = seq(0, 0.9, 0.05), pi0.method = "smoother", fdr.level = NULL, robust = FALSE, smooth.df = 3, smooth.log.pi0 = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
Vector of P values
}
  \item{lambda}{
 The value of the tuning parameter to estimate pi_0. Must be
          in [0,1). Optional.
}
  \item{pi0.method}{
 Either "smoother" or "bootstrap"; the method for
          automatically choosing tuning parameter in the estimation of
          pi_0, the proportion of true null hypotheses
}
  \item{fdr.level}{
 A level at which to control the FDR. Must be in (0,1].
          Optional; if this is selected, a vector of TRUE and FALSE is
          returned that specifies whether each q-value is less than
          fdr.level or not.}
  \item{robust}{
 An indicator of whether it is desired to make the estimate
          more robust for small p-values and a direct finite sample
          estimate of pFDR. Optional.
}
  \item{smooth.df}{
 Number of degrees-of-freedom to use when estimating pi_0
          with a smoother. Optional.
}
  \item{smooth.log.pi0}{
 If TRUE and 'pi0.method' = "smoother", pi_0 will be
          estimated by applying a smoother to a scatterplot of log pi_0
          estimates against the tuning parameter lambda. Optional.
}
  \item{\dots}{
Other potential arguments.
}
}

\author{
John D. Storey <jstorey@princeton.edu>
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p, lambda = seq(0, 0.9, 0.05), pi0.method = "smoother", 
    fdr.level = NULL, robust = FALSE, smooth.df = 3, smooth.log.pi0 = FALSE, 
    ...) 
{
    err.func <- "edge.qvalue"
    if (min(p) < 0 || max(p) > 1) {
        err.msg(err.func, "P-values not in valid range.")
        return(invisible(1))
    }
    if (length(lambda) > 1 && length(lambda) < 4) {
        err.msg(err.func, "If length of lambda greater than 1, you need at least 4 values.")
        return(invisible(1))
    }
    if (length(lambda) > 1 && (min(lambda) < 0 || max(lambda) >= 
        1)) {
        err.msg(err.func, "Lambda must be in [0,1).")
        return(invisible(1))
    }
    m <- length(p)
    if (length(lambda) == 1) {
        if (lambda < 0 || lambda >= 1) {
            err.msg(err.func, "Lambda must be in [0,1).")
            return(invisible(1))
        }
        pi0 <- mean(p >= lambda)/(1 - lambda)
        pi0 <- min(pi0, 1)
    }
    else {
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
            pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }
        if (pi0.method == "smoother") {
            if (smooth.log.pi0) {
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
                  pi0.boot[i] <- mean(p.boot > lambda[i])/(1 - 
                    lambda[i])
                }
                mse <- mse + (pi0.boot - minpi0)^2
            }
            pi0 <- min(pi0[mse == min(mse)])
            pi0 <- min(pi0, 1)
        }
        else {
            err.msg(err.func, "'pi0.method' must be one of 'smoother' or 'bootstrap'")
            return(invisible(1))
        }
    }
    if (pi0 <= 0) {
        err.msg(err.func, "The estimated pi0 <= 0. Check that you have valid\np-values or use another lambda method.")
        return(invisible(1))
    }
    if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
        1)) {
        err.msg(err.func, "'fdr.level' must be within (0,1].")
        return(invisible(1))
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
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, fdr.level = fdr.level, significant = (qvalue <= 
                fdr.level), lambda = lambda)
    }
    else {
        retval <- list(call = match.call(), pi0 = pi0, qvalues = qvalue, 
            pvalues = p, lambda = lambda)
    }
    class(retval) <- "qvalue"
    return(retval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
