\name{getSpanningSet}
\alias{getSpanningSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Split Data into Bins
}
\description{
Splits data into bins that span the range of intensities
}
\usage{
getSpanningSet(snm.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snm.obj}{
An object of class snm
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
function(snm.obj) {
  m <- rowMeans(snm.obj$M[snm.obj$nulls,])
  steps <- seq(quantile(snm.obj$M[snm.obj$nulls,],0.001),
               quantile(snm.obj$M[snm.obj$nulls,],0.999), length.out=snm.obj$nbins)
  th <- sapply(m, function(x) { 
    sum(x <= steps)
  }) 
  bins <- split(snm.obj$nulls, th)
  bins
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
