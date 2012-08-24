\name{calculate.nulls}
\alias{calculate.nulls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the Identity of the Null Probes
Identify Null Probes
}
\description{
Returns null probes given p -values and pi0 estimate.
}
\usage{
calculate.nulls(pvals, pi0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvals}{
Vector of P values
}
  \item{pi0}{
Pi0 of P values
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
function(pvals, pi0) {
  which(rank(pvals) > (1-pi0)*length(pvals))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
