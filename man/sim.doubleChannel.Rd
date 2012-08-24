\name{sim.doubleChannel}
\alias{sim.doubleChannel}
\title{
Simulated data for a double channel microarray experiment.
}
\description{
Simulates data used in snm function examples.
}
\usage{
sim.doubleChannel(seed)
}
\arguments{
  \item{seed}{
Numeric value used to seed random number generator.
}
}
\details{
Simulated data set influenced by a probe-specific biological variable, a probe-specific adjustment variables, and intensity-dependent
array and dye effects.  Data were simulated for a total of 25,000 probes and 20 arrays.  The biological variable is a dichotomous variable specifying two groups (Group 1 and Group 2), 
with 10 channels per group.  The continuous probe-specific adjustment variable is sampled from a Normal(1,0.1) and mimics an age effect.  
The baseline probe intensities  were sampled from a chi(1,2) distribution.  Any baseline intensities greater than 15 were set to a random value from the 
interval [15,16].  The random variation terms were sampled from a Normal(0,0.25) and the array and dye functions were defined by randomly sampling coefficients 
for a two-dimensional B-spline basis functions from a Normal(0,1).

Randomly selected subsets of 30\% and 20\% of the probes were defined as influenced by the biological groups and 
age variables, respectfully.  The magnitude of the biological effects were sampled from a Normal(1,0.3) distribution, and the probe-specific age effects 
from a Normal(1,0.1). An instance of this simulated data can be generated using the code in the examples section below. 
}
\value{
 	\item{\code{raw.data}}{a 25,000 by 50 matrix of simulated data generated according to the description above.}
  	\item{\code{true.nulls}}{a vector of indices corresponding to the rows in raw.data of the probes unaffected by the biological variable of interest} 
  	\item{\code{bio.var}}{a model matrix of the biological variable of interest.} 
  	\item{\code{adj.var}}{a model matrix of the adjustment variables}
  	\item{\code{int.var}}{a data frame of the intensity-dependent adjustment variables}
}
\author{
Brig Mecham <brig.mecham@sagebase.org>
}
\seealso{
\code{\link{snm}}, \code{\link{sim.singleChannel}}, \code{\link{sim.preProcessed}}, \code{\link{sim.refDesign}}
}
\examples{
\dontrun{
doubleChannel <- sim.doubleChannel(12346)
snm.obj <- snm(doubleChannel$raw.data, doubleChannel$bio.var, 
			doubleChannel$adj.var, doubleChannel$int.var)
ks.test(snm.obj$pval[doubleChannel$true.nulls],"punif")
}
}
\keyword{datagen}
