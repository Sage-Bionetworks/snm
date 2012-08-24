\name{construct.normalized.data}
\alias{construct.normalized.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Construct Normalized Data

}
\description{
This function removes all effects due to probe-specific adjustment variables.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
construct.normalized.data(snm.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snm.obj}{
An object of class snm.
}
}

\value{
Updated snm.obj.
}

\author{
Brig Mecham <brig.mecham@sagebase.org>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(snm.obj) {
  dat = snm.obj$dat
  n = snm.obj$n.arrays
  xx0 = snm.obj$adj.var
  df0 = snm.obj$df.null
  xx1 = cbind(snm.obj$adj.var, snm.obj$bio.var)
  P1 = solve(t(xx1) \%*\% xx1) \%*\% t(xx1)
  cfs1 = t(P1 \%*\% t(dat))
  res1 = dat - cfs1 \%*\% t(xx1)
  x1 <- cbind(rep(1,dim(dat)[2]), snm.obj$bio.var)
  cfs2 <- cfs1[, -(2:ncol(adj.var))]  
  fit1 <- cfs2 \%*\% t(x1)
  snm.obj$dat <- fit1 + res1
  snm.obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
