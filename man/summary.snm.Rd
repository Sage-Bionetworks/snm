\name{snm.summary}
\alias{snm.summary}
\alias{summary.snm}
\title{
Display summary information for an snm object
}
\description{
Provides a summary of the \code{snm} fit.
}
\usage{
snm.summary(object, cuts=c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1), \dots)
\method{summary}{snm}(object, \dots)
}
\arguments{
  \item{object}{
  An object of class "snm", i.e., output from the \code{snm} function.
  }
  \item{cuts}{
  Cut points at which to calculate the number of significant probes.
  }
  \item{\dots}{
  Optional arguments to send to the \code{\link[qvalue]{qvalue}} function in calculating statistical significance.
  }
}
\details{
A summary of the \code{snm} fit.

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
summary(snm.obj)
}
}
\keyword{misc}
