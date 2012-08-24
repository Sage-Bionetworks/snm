\name{buildBasisFunction}
\alias{buildBasisFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds Basis Fuction
}
\description{
Builds basis function used to evaluate RNA concentrations.
}
\usage{
buildBasisFunction(snm.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snm.obj}{
An object of class snm
}
}
\details{
Builds the basis spline function used to evaluate RNA concentrations.
}
\value{
Basis spline matrix.
}

\author{
Brig Mecham <brig.mecham@sagebase.org>
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (snm.obj) 
{
  lower <- min(rowMeans(snm.obj$dat))
  upper <- max(rowMeans(snm.obj$dat))
  basisSplineFunction <- ns(seq(lower,upper,length=100), df = snm.obj$spline.dim, 
        Boundary.knots = c(lower - 1, upper + 1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
