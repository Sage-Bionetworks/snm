\name{sim.singleChannel}
\alias{sim.singleChannel}
\title{
Simulate data from a single channel microarray experiment
}
\description{
Simulates single channel data used as an example for snm function call.
}
\usage{
sim.singleChannel(seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed}{Numeric value used to seed random number generator.}
}
\details{
Simulated data set influenced by a single probe-specific biological, two probe-specific adjustment variables, and intensity-dependent
array effects.  Data were simulated for a total of 25,000 probes and 50 arrays.  The biological variable is a dichotomous variable specifying two groups (Group 1 and Group 2), 
with 25 arrays sampled from each group.  The dichtomous probe-specific adjustment variables has 5 different levels and mimics a batch effect.  The 5 batches each contain
10 samples, and are balanced with respect to the biological grouping factor.  The continuous probe-specific adjustment variable is
sampled from a Normal(1,0.1) and mimics an age effect.  The baseline probe intensities  were sampled from a chi(1,2) distribution.  
Any baseline intensities greater than 15 were set to a random value from the 
interval [15,16].  The random variation terms were sampled from a Normal(0,0.25) and the array functions were defined by randomly sampling coefficients 
for a two-dimensional B-spline basis functions from a Normal(0,1).

Randomly selected subsets of 30\%, 10\%, and 20\% of the probes were defined as influenced by the biological groups, batch, and 
age variables, respectfully.  The magnitude of the biological effects were sampled from a Normal(1,0.3) distribution, the 
probe-specific batch effects from a Normal(0,0.3) and the probe-specific age effects from a Normal(1,0.1). An instance of 
this simulated data can be generated using the code in the examples section below. 

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
\code{\link{snm}}, \code{\link{sim.doubleChannel}}, \code{\link{sim.preProcessed}}, \code{\link{sim.refDesign}}
}
\examples{
\dontrun{
singleChannel <- sim.singleChannel(12345)
snm.obj <- snm(singleChannel$raw.data,
		      singleChannel$bio.var,
		      singleChannel$adj.var,
		      singleChannel$int.var)
ks.test(snm.obj$pval[singleChannel$true.nulls],"punif")
}
}
\keyword{datagen}
