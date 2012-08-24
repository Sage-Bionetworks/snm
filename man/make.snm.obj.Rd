\name{make.snm.obj}
\alias{make.snm.obj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Snm Object
}
\description{
Makes snm object and sets the values to those provided by the user. Does some quality control, but this should be expanded in the future.
}
\usage{
make.snm.obj(Y, bio.var, adj.var, int.var, spline.dim, nbins, weights,diagnose,rm.adj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Matrix of raw, log base 2 transformed data
}
  \item{bio.var}{
Model matrix of biological variables
}
  \item{adj.var}{
Model matrix of adjustment variables
}
  \item{int.var}{
Data frame with the different intensity-dependent variables presented along the columns
}
  \item{spline.dim}{
Dimension of basis spline used for modeling intensity-dependent effects
}
  \item{nbins}{
Number of bins to partition data into when estimating intensity-dependent effects
}
  \item{weights}{
Weights for each probe.  Values unchanged by the algorithm.  Used to down-weight influence of data when calculating intensity-dependent effects.
}
  \item{diagnose}{
A flag that determines if diagnostic plots should be produced.
}  
  \item{rm.adj}{
If set to TRUE, all adjustment variables are removed from the normalized data.  Set to FALSE if you intend to perform classical hypothesis on the normalized data.
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

function(Y, bio.var, adj.var,int.var,spline.dim, nbins,weights,diagnose,fast.snm,retBio) {
  if(colnames(bio.var)[1] == "(Intercept)") { bio.var <- bio.var[,-1] }
  else{ stop("Model matrix for biological variables must contain an Intercept term in the first column with the title (Intercept).  Please try again!")}
  snm.obj = list()
  snm.obj$n.arrays = ncol(Y)
  snm.obj$n.probes = nrow(Y)
  snm.obj$fast.snm = fast.snm
  snm.obj$retBio = retBio
  if(fast.snm) {
    snm.obj$adj.var <- cbind(rep(1,snm.obj$n.arrays))
    snm.obj$bio.var <- cbind(adj.var, bio.var)
    snm.obj$bio.var <- snm.obj$bio.var[,-1]
    snm.obj$user.adj.var <- adj.var
    snm.obj$user.bio.var <- bio.var
    snm.obj$df.full = ifelse(is.null(ncol(snm.obj$bio.var)), 1, ncol(snm.obj$bio.var)) + ncol(snm.obj$adj.var)
    snm.obj$df.null = ncol(snm.obj$adj.var)    
  }else{
    snm.obj$bio.var = bio.var
    snm.obj$adj.var = adj.var
    snm.obj$df.full = ifelse(is.null(ncol(bio.var)), 1, ncol(bio.var)) + ncol(adj.var)
    snm.obj$df.null = ncol(adj.var)
  }
  snm.obj$int.var = int.var
  snm.obj$individuals = NULL
  snm.obj$spline.dim = spline.dim
  snm.obj$nbins = nbins
  snm.obj$diagnose = diagnose
  snm.obj$dat=Y
  snm.obj$weights=weights
  snm.obj$r.dat=Y
  class(snm.obj) = "edge"
  snm.obj <- checkArguments(snm.obj)
  return(snm.obj)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
