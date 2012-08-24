\name{calcSig}
\alias{calcSig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs inference on data.
}
\description{
Calculates p and q values for a given model fit.  Accepts both biological and adjustment variables.
}
\usage{
calcSig(Y = NULL, bio.var = NULL, adj.var = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
Matrix of data
}
  \item{bio.var}{
Model matrix of biological variables
}
  \item{adj.var}{
Model matrix of adjustment variables
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class EDGE.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brig Mecham <brig.mecham@sagebase.org> and John D. Storey <jstorey@princeton.edu>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

seed <- 1234
sim.d1 <- sim.preProcessed(seed=seed,0.5,0.3,0.1)
sig <- calcSig(sim.d1$raw.data, sim.d1$bio.var, sim.d1$adj.var)
ks.test(sig$pval[sim.d1$true.nulls],"punif")$p

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
