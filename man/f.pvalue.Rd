\name{f.pvalue}
\alias{f.pvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates P values by fitting a full and reduced model to each probe 
}
\description{
See title
}
\usage{
f.pvalue(dat, mod, mod0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
data matrix
}
  \item{mod}{
full model
}
  \item{mod0}{
reduced model
}
}
\value{
Some p values
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
function(dat,mod,mod0){
  # This is a function for performing
  # parametric f-tests on the data matrix
  # dat comparing the null model mod0
  # to the alternative model mod. 
  n <- dim(dat)[2]
  m <- dim(dat)[1]
  df1 <- dim(mod)[2]
  df0 <- dim(mod0)[2]
  p <- rep(0,m)
  Id <- diag(n)
  
  resid <- dat \%*\% (Id - mod \%*\% solve(t(mod) \%*\% mod) \%*\% t(mod))
  resid0 <- dat \%*\% (Id - mod0 \%*\% solve(t(mod0) \%*\% mod0) \%*\% t(mod0))
  
  rss1 <- resid^2 \%*\% rep(1,n)
  rss0 <- resid0^2 \%*\% rep(1,n)
  
  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  p <-  1-pf(fstats,df1=(df1-df0),df2=(n-df1))
  return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
