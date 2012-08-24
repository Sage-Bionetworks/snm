\name{buildBasisSplineMatrix}
\alias{buildBasisSplineMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Build Basis Spline Matrix
}
\description{
Evaluations splines at RNA concentrations and returns their values
}
\usage{
buildBasisSplineMatrix(M, basisSplineFunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{
Matrix of pooled estimated RNA concentrations. The pooling strategy is defined by the nbins parameter. The element in position i,j is equal to the average RNA concentration on array j of the probes' whose average concentration is in the region spanned by bin i.
}
  \item{basisSplineFunction}{
B-spline basis function defined by buildBasisFunction
}
}

\author{
Brig Mecham <brig.mecham@sagebase.org>
}


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (M, basisSplineFunction) 
{
    lnp <- dim(M)[1]
    bSM <- apply(M, 2, function(x) {
        predict(basisSplineFunction, x)
    })
    bSM.model <- sapply(1:(dim(basisSplineFunction)[2]), function(x) {
        pos <- (x - 1) * lnp + 1
        as.numeric(bSM[pos:(pos + lnp - 1), ])
    })
    colnames(bSM.model) <- paste("Bt", 1:dim(bSM.model)[2], sep = "")
    bSM.model
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
