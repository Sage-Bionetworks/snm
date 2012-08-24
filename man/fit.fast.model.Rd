\name{fit.fast.model}
\alias{fit.fast.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fast Version of SNM Algorithm

}
\description{
Fits the Study-Specific Model without a mixed effects model.  This function estimates the intensity-dependent effects after down-weighting all of the probe-specific variables.  The default fit.model call down-weights only the biological variables.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fit.fast.model(obs.fit, snm.obj, basisSplineFunction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs.fit}{
  List of estimated coefficient matrices and residuals from full and reduced models
}
  \item{snm.obj}{
An object of class snm
}
  \item{basisSplineFunction}{
Basis spline function
}
}
\value{
Updated snm.obj with intensity-dependent effects. 
}

\author{
Brig Mecham <brig.mecham@sagebase.org>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


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
    Y.pooled[i,] = apply(matrix(snm.obj$r.dat[as.vector(bins[[i]]),],
              ncol=ncol(snm.obj$dat)),2,
              weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
    M.pooled[i,] = apply(matrix(snm.obj$M[as.vector(bins[[i]]),],
              ncol=ncol(snm.obj$M)),2,
              weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
  }

  BB <- predict(basisSplineFunction,M.pooled[,1])
  X <- kronecker(contr.sum(length(unique(snm.obj$int.var[,1])))[snm.obj$int.var[,1],], BB)
  for(i in 2:dim(snm.obj$int.var)[2]) {
    X <- cbind(X,
               kronecker(contr.sum(length(unique(snm.obj$int.var[,i])))[snm.obj$int.var[,i],], BB))
  }
  
  wts <- sapply(bins,length) / 10; wts[wts > 1] <- 1
  cfs <- summary(lm(as.numeric(t(scale(t(Y.pooled),scale=FALSE))) ~ -1+X,weights=rep(wts,times=snm.obj$n.arrays)))$coef[,1]

  beta = vector("list", dim(snm.obj$int.var)[2])
  beta[[1]] = matrix(cfs[1:(snm.obj$spline.dim * (length(unique(snm.obj$int.var[,1])) - 1))], ncol = length(unique(snm.obj$int.var[,1])) - 1)
  beta[[1]] = cbind(beta[[1]], -drop(beta[[1]] \%*\% rep(1, length(unique(snm.obj$int.var[,1])) - 1)))
  
  for(i in 2:(dim(snm.obj$int.var)[2])) {
    beta[[i]] = matrix(cfs[1:(snm.obj$spline.dim * (length(unique(snm.obj$int.var[,i])) - 1)) + snm.obj$spline.dim * (length(unique(snm.obj$int.var[,i-1]))- (i - 1))], 
      ncol = length(unique(snm.obj$int.var[,i])) - 1)
    beta[[i]] = cbind(beta[[i]], -drop(beta[[i]] \%*\% rep(1, length(unique(snm.obj$int.var[,i])) - 1)))
  }
  
  sapply(1:snm.obj$n.arrays, function(id) {
    preds <- predict(basisSplineFunction, snm.obj$M[,id])
    int.fx <- -preds \%*\% beta[[1]][,as.numeric(snm.obj$int.var[,1])[id]]
    for(i in 2:dim(snm.obj$int.var)[2]) {
      int.fx <- int.fx + -preds \%*\% beta[[i]][,as.numeric(snm.obj$int.var[,i])[id]]
    }
    -int.fx
  }) -> snm.obj$array.fx

  
# Add useful variables to snm.obj
  snm.obj$Y.pooled <- Y.pooled
  snm.obj$M.pooled <- M.pooled
  snm.obj$bin.densities <- sapply(bins,length)
  return(snm.obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
