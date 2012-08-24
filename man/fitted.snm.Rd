\name{snm.fitted}
\alias{snm.fitted}
\alias{fitted.snm}
\title{
Extract fitted values from an snm object
}
\description{
Computes fitted values under models used in \code{snm} normalization.
}
\usage{
snm.fitted(object, \dots)
\method{fitted}{snm}(object, \dots)
}
\arguments{
  \item{object}{
  Output from the \code{snm} function.
  }
  \item{\dots}{
  Not used.
  }
}
\details{
Returns the fitted values under the "null model" (adjustment variables only) and the "full model" (adjustment variables + biological variables). 
}
\value{
  \item{fit0}{
  Linear model fits when regression each probe's normalized data on the null model, \code{~adj.var}.
  }
  \item{fit1}{
  Linear model fits when regression each probe's normalized data on the full model, \code{~adj.var+bio.var}.  
  }
}
\references{
Mecham BH, Nelson PS, Storey JD (2010) Supervised normalization of
microarrays. Bioinformatics, 26: 1308-1315.
}
\author{
John D. Storey <jstorey@princeton.edu>
}
\note{
These fits are useful for investigating the quality of the study-specific model used in the normalization.  For example, the residuals can be obtained from the full model fit and examined for latent structure.}
\seealso{
\code{\link{snm}}, \code{\link{sim.singleChannel}}
}
\examples{
\dontrun{
singleChannel <- sim.singleChannel(12345)
snm.obj <- snm(singleChannel$raw.data,
		      singleChannel$bio.var,
		      singleChannel$adj.var[,-6],
		      singleChannel$int.var, num.iter=10)
snm.fit = fitted(snm.obj)
res1 = snm.obj$norm.dat - snm.fit$fit1
snm.svd = fast.svd(res1)
cor(snm.svd$v[,1], singleChannel$adj.var[,6])
plot(snm.svd$v[,1], singleChannel$adj.var[,6])
}
}
\keyword{misc}
