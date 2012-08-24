\name{sva}
\alias{sva}
\title{Estimate surrogate variables with an iterative algorithm from gene expression and model data}
\description{
  Estimate surrogate variables are estimated using either the iteratively re-weighted
  surrogate variable analysis algorithm of Leek and Storey (2008) or the two-step
  algorithm of Leek and Storey (2007).
}
\usage{
sva(dat, bio.var, adj.var=NULL, n.sv=NULL, num.iter=NULL, diagnose=TRUE, verbose=TRUE) 
}
\arguments{
  \item{dat}{Either an m genes by n arrays matrix of expression data or an object of class edge obtained from a previous sva function call.}
  \item{bio.var}{
A model matrix (see \code{\link[stats]{model.matrix}}) or data frame with \eqn{n} rows of the biological variables. If NULL, then all probes are treated as "null" in the algorithm. 
}
  \item{adj.var}{
A model matrix (see \code{\link[stats]{model.matrix}}) or data frame with \eqn{n} rows of the probe-specific adjustment variables.  If NULL, a model with an intercept term is used.
}
  \item{n.sv}{ Rank of dependence kernel.  If equal to NULL (default) this value is estimated from the data.}
  \item{num.iter}{ The number of iterations of the algorithm to perform.}
  \item{diagnose}{A flag telling the software whether or not to produce diagnostic output in the form of consecutive plots. TRUE produces the plot.} 
  \item{verbose}{A flag telling the software whether or not to display a report after each iteration. TRUE produces the output.}
}
\details{
  Surrogate variable estimates are formed based on unpublished modifications of the algorithms originally published
  in Leek and Storey (2007,2008). Surrogate variables can be included
  in a significance analysis to reduce dependence and confounding. 
}

\value{
  An object of class edge with the following values:  CURRENTLY THIS IS WRONG.  THE OUTPUT NEEDS TO BE CLEANED UP A BIT.
  \item{sv}{A n by n.sv matrix where each column is a distinct surrogate 
  variable (the main quantity of interest)}
  \item{pprob.gam}{A vector with the posterior probability estimates that
  each row is affected by dependence.} 
  \item{pprob.b}{A vector with the posterior probabiliity estimates that each
  row is affected by the variables in mod, but not in mod0.}
  \item{n.sv}{The number of suggorate variables estimated. }		
}

\references{
   Leek JT and Storey JD. (2008) A general framework for multiple testing 
   dependence. Proceedings of the National Academy of Sciences, 105: 
   18718-18723. \url{http://www.biostat.jhsph.edu/~jleek/publications.html}

   Leek JT and Storey JD. (2007) Capturing heterogeneity in gene expression 
   studies by surrogate variable analysis. PLoS Genetics, 3: e161. 
   \url{http://www.biostat.jhsph.edu/~jleek/publications.html}
}

\author{Brig Mecham \email{brig.mecham@sagebase.org}, John Storey \email{jstorey@princeton.edu}}
\seealso{\code{\link{snm}},\code{\link{edge}}, \code{\link{qvalue}}}
\examples{
  \dontrun{

seed <- 1234 
sim.d1 <- sim.preProcessed(seed=seed,0.5,0.3,0.1)

# Update and fit model 
sva.obj <- sva(sim.d1$raw.data, sim.d1$bio.var, NULL, n.sv=5,num.iter=5,diagnose=TRUE)
ps <- f.pvalue(sim.d1$raw.dat, model.matrix(~-1+sim.d1$bio.var+sva.obj$svd[[5]]$v), model.matrix(~sva.obj$svd[[5]]$v))
ks.test(ps[sim.d1$true.nulls],"punif")$p

# Update model and fit again
sva.obj2 <- sva(sva.obj,num.iter=5)
ps <- f.pvalue(sim.d1$raw.dat, model.matrix(~-1+sim.d1$bio.var+sva.obj2$svd[[10]]$v), model.matrix(~sva.obj2$svd[[10]]$v))
ks.test(ps[sim.d1$true.nulls],"punif")$p

# Now include one of the adjustment variables and fit
sva.obj <- sva(sim.d1$raw.data, sim.d1$bio.var, NULL, n.sv=5,num.iter=5,diagnose=TRUE)
ps <- f.pvalue(sim.d1$raw.dat, model.matrix(~-1+sim.d1$bio.var+sim.d1$adj.var[,6] + sva.obj$svd[[5]]$v), model.matrix(~sim.d1$adj.var[,6] + sva.obj$svd[[5]]$v))
ks.test(ps[sim.d1$true.nulls],"punif")$p


 }
}
\keyword{misc}
