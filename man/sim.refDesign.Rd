\name{sim.refDesign}
\alias{sim.refDesign}
\title{
Simulates data from a two-color microarray experiment using a reference design.
}
\description{
Simulates a two channel experiment with a reference design. Used as an example for snm function call.
}
\usage{
sim.refDesign(seed)
}
\arguments{
  \item{seed}{
Numeric value used to seed random number generator.
}
}
\details{
Simulated data set influenced by a probe-specific biological variable, and intensity-dependent
array and dye effects.  Data were simulated assuming a uniform reference design for a total of 25,000 probes and 20 arrays, each consisting of two channels.  
The reference channel consists of a single reference RNA population. The experimental channel measures a 
dichotomous biological variable specifying two groups (Group 1 and Group 2), with 10 samples per group.  The baseline probe intensities  
were sampled from a chi(1,2) distribution.  Any baseline intensities greater than 15 were set to a random value from the 
interval [15,16].  The random variation terms were sampled from a Normal(0,0.25) and the array and dye functions were defined by randomly sampling coefficients 
for a two-dimensional B-spline basis functions from a Normal(0,1).

A randomly selected subset of 30\% of the probes was defined as influenced by the biological group variable.  The magnitude of the biological 
effects were sampled from a Normal(1,0.3) distribution. An instance of this simulated data can be generated using the code in the examples section below. 
}
\value{
 	\item{\code{raw.data}}{a 25,000 by 50 matrix of simulated data generated according to the description above.}
  	\item{\code{true.nulls}}{a vector of indices corresponding to the rows in raw.data of the probes unaffected by the biological variable of interest} 
  	\item{\code{bio.var}}{a model matrix of the biological variable of interest.} 
  	\item{\code{adj.var}}{a model matrix of the probe-specific adjustment variables.  In this case set to NULL.}
  	\item{\code{int.var}}{a data frame of the intensity-dependent adjustment variables}
}
\author{
Brig Mecham <brig.mecham@sagebase.org> and John D. Storey <jstorey@princeton.edu>
}
\seealso{
\code{\link{snm}}, \code{\link{sim.singleChannel}}, \code{\link{sim.doubleChannel}}, \code{\link{sim.preProcessed}}
}
\examples{
\dontrun{
# reference design on channel level data
# reference channel is included in bio.var
refChannel <- sim.refDesign(12347)
snm.obj <- snm(refChannel$raw.data,refChannel$bio.var, refChannel$adj.var, refChannel$int.var)
ks.test(snm.obj$pval[refChannel$true.nulls],"punif")

# reference design on log ratio data
refChannel <- sim.refDesign(12347)
refChannel$raw.data = refChannel$raw.data[,1:20]-refChannel$raw.data[,21:40]
# removing reference channel variables
refChannel$bio.var = refChannel$bio.var[1:20,-3]
refChannel$int.var = data.frame(refChannel$int.var[1:20,1])
snm.obj <- snm(refChannel$raw.data,refChannel$bio.var, refChannel$adj.var, refChannel$int.var)
ks.test(snm.obj$pval[refChannel$true.nulls],"punif")
}
}
\keyword{datagen}
