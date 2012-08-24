\name{get.resids}
\alias{get.resids}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate Residuals 
}
\description{
Get Residuals from data given a model
}
\usage{
get.resids(Y, mod)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Matrix of data
}
  \item{mod}{
Model matrix parameterizing coefficients
}
}



\author{
Brig Mecham <brig.mecham@sagebase.org> and John D. Storey <jstorey@princeton.edu>
}




\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data.matrix, design.matrix) 
{
    x.qr <- qr(design.matrix)
    x.det <- prod(diag(x.qr$qr))
    x.cov.unscaled <- solve(t(qr.R(x.qr)) \%*\% qr.R(x.qr))
    x.df <- nrow(design.matrix) - x.qr$rank
    Fit.coef <- qr.coef(x.qr, t(data.matrix))
    Fit.res <- qr.resid(x.qr, t(data.matrix))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
