\name{sim.function.var}
\alias{sim.function.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function simulates data as a function of a variable
}
\description{
This function simulates data that is a function of a variable.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim.function.var(data = NULL, variable = NULL, rows = NULL, sample.from = NULL, spline.dim = NULL)
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
  \item{spline.dim}{
Dimension of spline for function to be simulated.
}
}
\value{
Matrix of dimension equal to that of data containing effects of the simulated function.  
}

\author{
Brig Mecham <brig.mecham@sagebase.org>
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(data=NULL,variable=NULL,rows=NULL, sample.from=NULL, spline.dim=NULL) {

#  Make certain class of variable is numeric
  if(class(variable) != "numeric" & class(variable) != "integer") { 
    stop("Class of variable must be numeric.")
  }
  
  Z <- ns(variable, df=spline.dim)
# Handle the rows variable.  If its a single number, then its a proportion.  Otherwise, its a vector of probe indices to be modified.
  lab <- rep(0, nrow(data))  # Create vector of dummy variables    
  if(length(rows)==1) {
    if(rows > 1) { stop("Proportion of data to be influenced by variable greater than 1.")}
    these <- sample(1:nrow(data), nrow(data) * rows) # Sample rows to be modified
    lab[these] <- 1  #Identify rows to be modified
  }else{
    if(max(rows) > nrow(data)) { stop("At least one row of data matrix to be influenced by variable is larger than total number of rows in data")}
    if(min(rows) < 0) { stop("Rows to be influenced by variable must be greater than 0")}
    lab[rows] <- 1  # Identify rows to be modified
  }

  x <- model.matrix(~-1+Z)  # Create model matrix for variable
  sample.this.many <- sum(lab==1) * ncol(x)  # Count how many probes to be modified
  cfs <- matrix(0, nr=length(lab), nc=ncol(x))  #Initialize matrix of coefficients

  if(is.list(sample.from)) {
    for(spd in 1:length(sample.from)) { 
                                        # Estimate coefficients and add to coefficients matrix
      cfs[which(lab==1),spd] <- do.call(sample.from[[spd]]$func, as.list(c(n=sum(lab==1), sample.from[[spd]]$params)))
    }
  }else{
    stop("Spline coefficients must be a list")
  }

  cfs.mat <- cfs \%*\% t(x)  #  Estimate overall effect
  t(apply(cfs.mat,1,function(x) {  #Make all positive.  Don't want streaks in scatter plots
    if(min(x) < 0) {
      x - min(x)
    }else{
      x
    }
  })) -> cfs.mat
  cfs.mat  # Return effects
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
