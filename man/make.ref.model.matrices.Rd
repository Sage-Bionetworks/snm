\name{make.ref.model.matrices}
\alias{make.ref.model.matrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Model Matrices for Random Effects Model
}
\description{
Makes model matrices for random effects models
}
\usage{
make.ref.model.matrices(snm.obj, exp)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snm.obj}{
An object of class snm
}
  \item{exp}{
New environment to hold the lmer formula
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
function (snm.obj, exp) 
{
    F.mats <- list()
    spline.re.model <- ""
    for (i in 1:length(snm.obj$int.var)) {
        ide <- names(snm.obj$int.var)[i]
        spline.re.model <- paste(spline.re.model, "+ ( Bt", 1, 
            "-1 | ", ide, ")", sep = "")
        for (d in 2:snm.obj$spline.dim) {
            spline.re.model <- paste(spline.re.model, "+ ( Bt", 
                d, "-1 | ", ide, ")", sep = "")
        }
        F.mats[[i]] <- model.matrix(~-1 + snm.obj$int.var[, 
            i])
    }
    if (dim(snm.obj$adj.var)[2] > 0) {
        ref.model <- paste( "probes+", paste("probes:", colnames(snm.obj$adj.var)[-1],
            collapse = "+"),sep="+" )
    }
    else {
        ref.model <- "probes + "
    }
    ref.model <- as.formula(paste("y ~ -1 +", ref.model, spline.re.model), 
        env = exp)
    names(F.mats) <- names(snm.obj$int.var)
    list(ZF = ref.model, F.mats = F.mats)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
