\name{snm}
\alias{snm}
\title{
Perform a supervised normalization of microarray data
}
\description{
Implement Supervised Normalization of Microarrays on a gene
  expression matrix.  Requires a set of biological covariates of
  interest and at least one probe-specific or intensity-dependent 
  adjustment variable.
}
\usage{
snm(raw.dat, bio.var=NULL, adj.var=NULL, int.var=NULL,
    weights=NULL, spline.dim = 4, num.iter = 10, nbins=50,
    rm.adj=FALSE, verbose=TRUE, diagnose=TRUE)
}
\arguments{
  \item{raw.dat}{
An \eqn{m} probes by \eqn{n} arrays matrix of expression data.  If the user wishes to remove intensity-dependent effects, then we request
the matrix corresponds to the raw, log transformed data.
}
  \item{bio.var}{
A model matrix (see \code{\link[stats]{model.matrix}}) or data frame with \eqn{n} rows of the biological variables. If NULL, then all probes are treated as "null" in the algorithm. 
}
  \item{adj.var}{
A model matrix (see \code{\link[stats]{model.matrix}}) or data frame with \eqn{n} rows of the probe-specific adjustment variables.  If NULL, a model with an intercept term is used.
}
  \item{int.var}{
A data frame with \eqn{n} rows of type factor with the unique levels of intensity-dependent effects. Each column parametrizes a unique source of intensity-dependent effect (e.g., array effects for column 1 and dye effects for column 2).
}
  \item{weights}{
A vector of length \eqn{m}.  Values unchanged by algorithm, used to control the influence of each probe on the intensity-dependent array effects.
}
  \item{spline.dim}{
Dimension of basis spline used for array effects.
}
  \item{num.iter}{
Number of iterations to run.
}
  \item{nbins}{
Number of bins used by binning strategy.  Array effects are calculated from a \eqn{nbins} x \eqn{n} data matrix, where the \eqn{(i,j)} value is equal to that bin \eqn{i}'s average intensity on array \eqn{j}.
}
  \item{rm.adj}{
If set to FALSE, then only the intensity dependent effects have been removed from the normalized data, implying the effects from the adjustment variables are still present.  If TRUE, then the adjustment variables effects and the intensity dependent effects are both removed from the returned normalized data.
}
  \item{verbose}{
A flag telling the software whether or not to display a report after each iteration. TRUE produces the output.
}
  \item{diagnose}{
A flag telling the software whether or not to produce diagnostic output in the form of consecutive plots. TRUE produces the plot.
} 
}
\details{
This function implements the supervised normalization of microarrays algorithm described in Mecham, Nelson, and Storey (2010).
}
\value{
  \item{norm.dat}{
The matrix of normalized data.  The default setting is rm.adj=FALSE, which means that only the intensity-dependent effects have been subtracted from the data.  If the user wants the adjustment variable effects removed as well, then set rm.adj=TRUE when calling the \code{snm} function.
}
  \item{pvalues}{
A vector of p-values testing the association of the biological variables with each probe.  These p-values are obtained
from an ANOVA comparing models where the full model contains both the probe-specific biological and adjustment variables versus a reduced
model that just contains the probe-specific adjustment variables.  The data used for this comparison has the intensity-dependent variables
removed.  These returned p-values are calculated after the final iteration of the algorithm.
}
  \item{pi0}{
The estimated proportion of true null probes \eqn{pi_0}, calculated after the final iteration of the algorithm.
}
  \item{iter.pi0s}{
A vector of length equal to num.iter containing the estimated \eqn{pi_0} values at each iteration of the \code{snm} algorithm.  These values should converge and any non-convergence suggests a problem with the data, the assumed model, or both
}
  \item{nulls}{
A vector indexing the probes utilized in estimating the intensity-dependent effects on the final iteration.
}
  \item{M}{
A matrix containing the estimated probe intensities for each array utilized in estimating the intensity-dependent effects on the final iteration.  For memory parsimony, only a subset of values spanning the range is returned, currently nbins*100 values.
}
  \item{array.fx}{
A matrix of the final estimated intensity-dependent array effects.  For memory parsimony, only a subset of values spanning the range is returned, currently nbins*100 values.
}
  \item{bio.var}{
The processed version of the same input variable.  
}
  \item{adj.var}{
The processed version of the same input variable.  
}
  \item{int.var}{
The processed version of the same input variable.  
}
  \item{df0}{
Degrees of freedom of the adjustment variables.
}
  \item{df1}{
Degrees of freedom of the full model matrix, which includes the biological variables and the adjustment variables.
}
  \item{raw.dat}{
The input data.
}
  \item{rm.var}{
Same as the input (useful for later analyses).
}
  \item{call}{
Function call.
}
}
\references{
Mecham BH, Nelson PS, Storey JD (2010) Supervised normalization of
microarrays. Bioinformatics, 26: 1308-1315.
}
\author{
Brig Mecham <brig.mecham@sagebase.org> and John D. Storey <jstorey@princeton.edu>
}
\note{
It is necessary for adj.var and adj.var+bio.var to be valid model matrices (e.g., the models cannot be over-determined).  

We suggest that the probe level data be analyzed on the log-transformed
scale, particularly if the user wishes to remove intensity-dependent
effects. It is recommended that the normalized data (and resulting
inference) be inspected for latent structure using Surrogate Variable
Analysis (Leek and Storey 2007, PLoS Genetics).
}
\seealso{
\code{\link[stats]{model.matrix}}, \code{\link{plot.snm}}, \code{\link{fitted.snm}}, \code{\link{summary.snm}}, \code{\link{sim.singleChannel}}, \code{\link{sim.doubleChannel}}, \code{\link{sim.preProcessed}}, \code{\link{sim.refDesign}}
}
\examples{
singleChannel <- sim.singleChannel(12345)
snm.obj <- snm(singleChannel$raw.data,
		      singleChannel$bio.var,
		      singleChannel$adj.var,
		      singleChannel$int.var)
ks.test(snm.obj$pval[singleChannel$true.nulls],"punif")
plot(snm.obj)
summary(snm.obj)
snm.fit = fitted(snm.obj)
}
\keyword{htest}
\keyword{models}
