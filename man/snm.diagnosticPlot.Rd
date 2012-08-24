\name{snm.diagnostic.plot}
\alias{snm.diagnostic.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Produce Diagnostic Plots For SNM Algorithm
}
\description{
This function produces plots useful for diagnosing problems with an assumed study-specific model.
}
\usage{
snm.diagnostic.plot(obs.fit, obs.stat, snm.obj, basisSplineFunction, n.iter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs.fit}{
  List of estimated coefficient matrices and residuals from full and reduced models
  }
  \item{obs.stat}{
  List of statistic obtained from model fit.  
}
  \item{snm.obj}{
  An object of class snm.
}
  \item{basisSplineFunction}{
  Basis spline functions for intenisty-dependent functions
}
  \item{n.iter}{
  Current iteration number
  }
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
function (obs.fit,obs.stat,snm.obj,
          basisSplineFunction)
{
  par(mfrow=c(2,2))
  bds <- ifelse(snm.obj$bin.densities > 500, 500, snm.obj$bin.densities)
  plot(rowMeans(snm.obj$M.pooled),
       bds,ylim=c(0,500),
       xlab="Estimated RNA Concentration", ylab="Probes per bin", main="Null probes per bin")

  u <- fast.svd(obs.fit$res1,tol=0);
  plot(u$d^2 / sum(u$d^2), ylim=c(0,1),
       main="Latent Structure",
       xlab="Surrogate Variable",
       ylab="Percent Variance Explained")

  if(dim(snm.obj$dat)[1] > 5000) {
    th <- sample(dim(snm.obj$dat)[1], 5000)
  }else{
    th <- 1:dim(snm.obj$dat)[1]
  }

  oo <- order(snm.obj$M[th,1])
  plot(snm.obj$M[th[oo],1], snm.obj$array.fx[th[oo],1],
       type="l",lwd=3, ylim=range(snm.obj$array.fx[th[oo],]),
       xlab="Estimated RNA Concentration", ylab="Array Effect",main="Estimated Array Effects")
  sapply(2:dim(snm.obj$r.dat)[2],function(id) {
    oo <- order(snm.obj$M[th,id])
    points(snm.obj$M[th[oo],id], snm.obj$array.fx[th[oo],id],type="l",lwd=3)
  }) -> hmm

  hist(obs.stat$pval, xlab="P values",main="P value distribution");
  abline(v=sort(obs.stat$pval)[length(obs.stat$pval)*(1-snm.obj$pi0)],col="red")

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
