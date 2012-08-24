\name{snm.plot}
\alias{snm.plot}
\alias{plot.snm}
\title{
Display plots for an snm object
}
\description{
Creates a diagnostic plot of the \code{snm} fit.
}
\usage{
snm.plot(x, col.by=NULL, \dots)
\method{plot}{snm}(x, \dots)
}
\arguments{
  \item{x}{
  Output from the \code{snm} function.
  }
  \item{col.by}{
  A factor vector of length equal to the number of arrays providing a grouping by which to color the intensity-dependent effects.  Instead of a factor vector, the input may also be a model matrix composed only of 0's and 1's with the number of rows equal to the number of arrays and number of columns less than or equal to the number of arrays.
  }
  \item{...}{
  Arguments passed to the plot functions.  Not recommended.
  }
}
\details{
A four panel plot composed of the following:
\enumerate{
\item Convergence of \eqn{\pi_0}{pi_0} estimates over model fitting iterations.  The \eqn{pi_0}{pi_0} estimates for each iteration are compared to the \eqn{pi_0} estimate calculated during the final model fit.
\item A scree plot of the principal components analysis of the full model residual matrix. 
\item A plot of the estimated intensity-dependent effects. 
\item A histogram of the p-values testing each probe for an asssociation with the biological variables (\code{bio.var}). All probes to the right of the vertical red line are the least \eqn{\hat{\pi}_0}{pi0.hat} significant probes (i.e., those used in estimating intensity-dependent effects).  The dashed horizontal line is the \eqn{pi_0}{pi_0} estimate from the final model fit. 
}
}
\value{
Nothing of interest.
}
\references{
Mecham BH, Nelson PS, Storey JD (2010) Supervised normalization of
microarrays. Bioinformatics, 26: 1308-1315.
}
\author{
John D. Storey <jstorey@princeton.edu>
}
\seealso{
\code{\link{snm}}, \code{\link{sim.singleChannel}}
}
\examples{
\dontrun{
singleChannel <- sim.singleChannel(12345)
snm.obj <- snm(singleChannel$raw.data,
		      singleChannel$bio.var,
		      singleChannel$adj.var,
		      singleChannel$int.var, num.iter=10)
plot(snm.obj, col.by=snm.obj$bio.var) #color by biological group
plot(snm.obj, col.by=snm.obj$adj.var[,-6]) #color by batch
}
}
\keyword{misc}
