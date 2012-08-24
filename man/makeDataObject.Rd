\name{makeDataObject}
\alias{makeDataObject}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Makes data object for fitting the random effects model
}
\description{
Makes data object for fitting the random effects model
}
\usage{
makeDataObject(Y.pooled, np, snm.obj, exp, bins)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y.pooled}{
Matrix of pooled raw data. Element i,j is the average raw data for every probe in bin i on array j.
}
  \item{np}{
list matching bins to probes
}
  \item{snm.obj}{
An object of class snm
}
  \item{exp}{
New environment to hold the lmer formula
}
  \item{bins}{
list matching bins to probes
}
}



\author{
Brig Mecham <brig.mecham@sagebase.org> and John D. Storey <jstorey@princeton.edu>
}




\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y.pooled, np, snm.obj, bSM.model, exp,bins) 
{
   lnp <- length(np)
   num.arrays <- dim(Y.pooled)[2]
   signal <- data.frame(y = as.numeric(Y.pooled), probes = factor(rep(1:lnp, 
                                             times = num.arrays)))
   f <- lapply(1:length(snm.obj$int.var), function(x) {
     i <- snm.obj$int.var[, x]
     rep(i, each = lnp)
   })
   f <- as.data.frame(f)
   colnames(f) <- names(snm.obj$int.var)
   if (dim(snm.obj$adj.var)[2] > 1) {
     z <- lapply(2:(dim(snm.obj$adj.var)[2]), function(x) {
       i <- snm.obj$adj.var[, x]
       rep(i, each = lnp)
     })
     z <- as.data.frame(z)
     colnames(z) <- colnames(snm.obj$adj.var)[-1]
     df <- cbind(cbind(cbind(signal, z), f), bSM.model)
   }else {
     df <- cbind(cbind(signal, f), bSM.model)
   }
   df$weights = sapply(bins,length) / 50 #mean(sapply(bins,length))
   df$weights[df$weights >1] <- 1
   TMP <- sapply(1:dim(df)[2], function(x) {
     assign(colnames(df)[x], df[, x], envir = exp)
   })
   as.data.frame(df)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
