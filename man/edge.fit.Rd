\name{edge.fit}
\alias{edge.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fits Model to Data
}
\description{
Fits model to data and returns many objects useful when performing hypotheses tests.
}
\usage{
edge.fit(edge.obj, odp = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edge.obj}{
Object of class edge
}
  \item{odp}{
Flag determining whether or not to implement ODP on P-values.
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
function (edge.obj, odp = FALSE) 
{
    err.func = "edge.fit"
    dat = edge.obj$dat
    n = edge.obj$n.arrays
    xx0 = edge.obj$adj.var
    df0 = edge.obj$df.null
    if (!odp) {
        xx1 = cbind(edge.obj$adj.var, edge.obj$bio.var)
    }
    else {
        xx1 = edge.obj$bio.var
    }
    df1 = edge.obj$df.full
    if (!is.null(edge.obj$ind)) {
        xxi = edge.obj$ind
        Hi = xxi \%*\% solve(t(xxi) \%*\% xxi) \%*\% t(xxi)
        fit.ind = t(Hi \%*\% t(dat))
        dat = dat - fit.ind
        xx1 = xx1 - Hi \%*\% xx1
        xx0 = xx0 - Hi \%*\% xx0
        xx1 = rm.zero.cols(xx1)
        xx0 = rm.zero.cols(xx0)
    }
    H0 = xx0 \%*\% solve(t(xx0) \%*\% xx0) \%*\% t(xx0)
    fit0 = t(H0 \%*\% t(dat))
    res0 = dat - fit0
    if (odp) {
        xx1 = xx1 - H0 \%*\% xx1
        H1 = xx1 \%*\% solve(t(xx1) \%*\% xx1) \%*\% t(xx1)
        fit1 = t(H1 \%*\% t(res0))
        res1 = res0 - fit1
        var1 = rowSums(res1^2)/(n - df1)
        var0 = rowSums(res0^2)/(n - df0)
        return(list(fit1 = fit1, fit0 = fit0, res1 = res1, res0 = res0, 
            var1 = var1, var0 = var0, ind = edge.obj$ind, odp = odp))
    }
    else {
        H1 = xx1 \%*\% solve(t(xx1) \%*\% xx1) \%*\% t(xx1)
        fit1 = t(H1 \%*\% t(dat))
        res1 = dat - fit1
        return(list(fit0 = fit0, res1 = res1, res0 = res0, ind = edge.obj$ind, 
            odp = odp))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
