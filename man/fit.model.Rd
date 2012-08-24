\name{fit.model}
\alias{fit.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds and Fits the Mixed Effects Model
Fit Random Effects Model
}
\description{
Function estimates the binned matrices and fits the random effects models
Function designed to fit the random effects model to the data.  
}
\usage{
fit.model(obs.fit, snm.obj, basisSplineFunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{obs.fit}{
List of estimated coefficient matrices and residuals from full and reduced models
}
\item{snm.obj}{ An object of class snm}

\item{basisSplineFunction}{Basis spline function}
}



\author{
Brig Mecham <brig.mecham@sagebase.org> and John D. Storey <jstorey@princeton.edu>
}




\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obs.fit,snm.obj,
          basisSplineFunction)
{
  snm.obj$M <- snm.obj$dat - obs.fit$res1
  snm.obj$M[snm.obj$nulls,] <- obs.fit$fit0[snm.obj$nulls,]

# Split the data into nbins bins based on their mean intensities
  bins <- getSpanningSet(snm.obj)
  
# Build the matrix of weighted raw data and matrix of weighted fitted values for each bin.
  lnp <- length(bins)
  np <- 1:lnp
  Y.pooled <- 0*snm.obj$dat[np,]
  M.pooled <- 0*snm.obj$M[np,]
  for(i in 1:lnp) {
    Y.pooled[i,] = apply(matrix(snm.obj$r.dat[as.vector(bins[[i]]),], ncol=ncol(snm.obj$dat)),2,
              weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
    M.pooled[i,] = apply(matrix(snm.obj$M[as.vector(bins[[i]]),], ncol=ncol(snm.obj$M)),2,
              weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
  }
  
# Build the basis spline matrix for the pooled coefficients.
  bSM.model <- buildBasisSplineMatrix(M.pooled, basisSplineFunction)
  exp <- new.env()
# Build the data object and fit the mixed effects model
  expObj <- makeDataObject(Y.pooled, np, snm.obj,
                           bSM.model, exp,bins)
  model.objects <- make.ref.model.matrices(snm.obj, exp)
  rff <- do.call("lmer", list(model.objects$ZF,expObj,NULL,TRUE,list(),NULL,FALSE,TRUE,TRUE,expObj$weights))
# Add useful variables to snm.obj
  snm.obj$E.pooled <- matrix(rff@resid, nr=dim(Y.pooled)[1])
  snm.obj$Y.pooled <- Y.pooled
  snm.obj$M.pooled <- M.pooled
  snm.obj$bin.densities <- sapply(bins,length)
  qs <- list(pi0=snm.obj$pi0, pvalues=snm.obj$pvalues)
#  Parse the lmer object!!! 
  snm.obj$array.fx <- calcArrayEffects(ranef(rff), basisSplineFunction, 
                               snm.obj,
                               model.objects,snm.obj$M)
  return(snm.obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
