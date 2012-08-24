\name{rm.zero.cols}
\alias{rm.zero.cols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cleans Model Matrices
}
\description{
Removes Zero Columns
}
\usage{
rm.zero.cols(xx, eps = 10^(-12))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{
Model matrix
}
  \item{eps}{
Cutoff
}
}
\author{
John D. Storey <jstorey@princeton.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xx, eps = 10^(-12)) 
{
    n = ncol(xx)
    vv = NULL
    for (i in 1:n) {
        if (sum(abs(xx[, i])) < eps) {
            vv = c(vv, i)
        }
    }
    if (is.null(vv)) {
        return(xx)
    }
    else {
        return(xx[, -vv])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
