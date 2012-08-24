\name{sim.intensity.dep}
\alias{sim.intensity.dep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Intensity Dependent Variables

}
\description{
This function simulates intensity dependent variables.  
}
\usage{
sim.intensity.dep(data = NULL, variable = NULL, spline.dim = NULL, sample.from = NULL, rows = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
The data matrix of simulated data the user wishes to adjust.
}
  \item{variable}{
The variable the user wishes to introduce.
}
  \item{spline.dim}{
Dimension of intensity dependent function
}
  \item{sample.from}{
A list specifying a function that coefficients can be sampled from (e.g. rnorm, rchisq, etc.) and the relevant parameters (e.g. mean, sd, df, etc).  The first item in the list should be named func and the second named params.
}
  \item{rows}{
Either a number between 0 and 1 or a vector of row IDs.  A number between 0 and 1 specifies the proportion of data to be influenced by the variable.  A vector specifies the rows to be influenced.   Equal to 1 by default.
}
}
\value{
Matrix of dimension equal to that of data containing the intensity-dependent effects.  
}

\author{
Brig Mecham <brig.mecham@sagebase.org>
}


\examples{
data <- matrix(rchisq(100,1,2),nr=1000,nc=12)
int.var <- data.frame(array=factor(1:12))
array.effect <- sim.intensity.dep(data, int.var$array, 2, list(func=rnorm, params=c(mean=0,sd=0.5)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
