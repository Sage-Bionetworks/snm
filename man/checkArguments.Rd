\name{checkArguments}
\alias{checkArguments}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Not very useful quite yet, but if software is expanded to include other types of data this function will become more relevant.
Verify User Input
}
\description{
 Not very useful quite yet, but if software is expanded to include other types of data this function will become more relevant.
At this time a rather useless function.  Purpose is mainly for any extensions of this package that might require more complex input (e.g. CNV arrays, surface modeling, etc).
}
\usage{
checkArguments(edge.obj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edge.obj}{
An object of class snm.
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
function (edge.obj) 
{
    m <- dim(edge.obj$dat)[1]
    n <- dim(edge.obj$dat)[2]
    d <- dim(edge.obj$bio.var)[2]
    n.bc <- dim(edge.obj$bio.var)[1]
#    if (n.bc != n) {
#        print("Number of columns of Y and number of rows of bio.vars are not equal\n ")
#        stop()
#    }
#    n.r.c <- dim(edge.obj$adj.var)[1]
#    if (n.r.c != n && dim(edge.obj$adj.var)[1] != 0) {
#        print("Number of columns of Y and number of rows of technical.covariates are not equal\n ")
#        stop()
#    }
#    if (!is.matrix(edge.obj$dat)) {
#        print("Y is not a numeric matrix\n")
#        stop()
#    }
    return(edge.obj)    
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
