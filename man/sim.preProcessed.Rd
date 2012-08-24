\name{sim.preProcessed}
\alias{sim.preProcessed}
\title{
Simulate data from a microarray experiment without any intensity-dependent effects.
}
\description{
Simulated data set influenced by a single probe-specific biological and two probe-specific adjustment variables.  This parameters for this data
are identical to single channel simulated data available as sim.singleChannel(seed) with the difference that this example does not include
the intensity-dependent effects. Consult the corresponding help file for details on this simulation.
}
\usage{
sim.preProcessed(seed,perc.bio=0.3,perc.batches=0.3,perc.height=0.1,np=50000)
}
\arguments{
  \item{seed}{Numeric value used to seed random number generator.}
  \item{perc.bio}{Percentage of probes influenced by biological variable.}
  \item{perc.batches}{Percentage of probes influenced by batch.}
  \item{perc.height}{Percentage of probes influenced by height.}
  \item{np}{Number of probes to simulate.}
}

\value{
    \item{\code{raw.data}}{a 25,000 by 50 matrix of simulated data generated according to the description above.}
  	\item{\code{true.nulls}}{a vector of indices corresponding to the rows in raw.data of the probes unaffected by the biological variable of interest} 
  	\item{\code{bio.var}}{a model matrix of the biological variable of interest.} 
  	\item{\code{adj.var}}{a model matrix of the adjustment variables}
  	\item{\code{int.var}}{set to NULL}
}
\author{
Brig Mecham <brig.mecham@sagebase.org> 
}
\seealso{
\code{\link{snm}}, \code{\link{sim.doubleChannel}}, \code{\link{sim.singleChannel}}, \code{\link{sim.refDesign}}
}
\examples{
preProcessed <- sim.preProcessed(12347)
snm.obj <- snm(preProcessed$raw.data, 
                      preProcessed$bio.var,
                      preProcessed$adj.var, rm.adj=TRUE)
ks.test(snm.obj$pval[preProcessed$true.nulls],"punif")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
