\name{edge.glr}
\alias{edge.glr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Test Statistics and P values
}
\description{
Calculates test statistics and P values from obs.fit, which is returned from edge.fit. 
}
\usage{
edge.glr(obs.fit, df1, df0, norm.pval = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs.fit}{
Object returned from edge.fit
}
  \item{df1}{
Degrees of freedom of full model
}
  \item{df0}{
Degrees of freedom of reduced model
}
  \item{norm.pval}{
Flage determining whether or not P value should be normalized.
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
function (obs.fit, df1, df0, norm.pval = FALSE) 
{
    err.func = "edge.glr"
    res1 = obs.fit$res1
    res0 = obs.fit$res0
    n = ncol(res1)
    rss1 = apply(res1, 1, function(x) {
        sum(x^2)
    })
    rss0 = apply(res0, 1, function(x) {
        sum(x^2)
    })
    stat = ((rss0 - rss1)/(df1 - df0))/(rss1/(n - df1))
    if (norm.pval) {
        pval = 1 - pf(stat, df1 = (df1 - df0), df2 = (n - df1))
        return(list(stat = stat, pval = pval))
    }
    else {
        return(list(stat = stat, pval = NULL))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
