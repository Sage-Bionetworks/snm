\name{sim.probe.specific}
\alias{sim.probe.specific}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a Probe-Specific Variable
}
\description{
This function is designed to simulate a probe-specific variable.  It accepts one variable at a time.  The variable need be either a factor or numeric vectors, any character vector will be modeled as a factor.
}
\usage{
sim.probe.specific(data = NULL, variable = NULL, rows = NULL, sample.from = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data matrix of simulated data the user wishes to adjust.
}
  \item{variable}{
The variable the user wishes to introduce.
}
  \item{rows}{
Either a number between 0 and 1 or a vector of row IDs.  A number between 0 and 1 specifies the proportion of data to be influenced by the variable.  A vector specifies the rows to be influenced. 
}
  \item{sample.from}{
A list specifying a function that coefficients can be sampled from (e.g. rnorm, rchisq, etc.) and the relevant parameters (e.g. mean, sd, df, etc).  The first item in the list should be named func and the second named params.
}
}
\value{
Matrix of dimension equal to that of data containing the probe-specific effects.}

\author{
Brig Mecham <brig.mecham@sagebase.org>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

data <- matrix(rchisq(100,1,2),nr=1000,nc=12)
bio.var <- data.frame(groups=rep(c("A","B"),each=6))
group.effect <- sim.probe.specific(data, bio.var$groups, 0.3, list(func=rnorm,params=c(mean=1,sd=0.3)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ simulations }
\keyword{ microarray }% __ONLY ONE__ keyword per line
