\name{err.msg}
\alias{err.msg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Error Message
}
\description{
Reports and Error Message
}
\usage{
err.msg(err.func = "edge", msg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{err.func}{
Variable that declares the function in which the error occurs.
}
  \item{msg}{
Message to be reported.
}
}
\author{
John Storey <jstorey@princeton.edu>
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (err.func = "edge", msg) 
{
    cat("\n")
    cat("\t")
    cat("ERROR in the", err.func, "function: ", "\n")
    cat("\t", msg, "\n\n")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
