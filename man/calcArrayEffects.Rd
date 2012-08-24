\name{calcArrayEffects}
\alias{calcArrayEffects}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calulate Array Effects Using Random Effects Model
}
\description{
Fits the mixed effects model
}
\usage{
calcArrayEffects(rff, basisSplineFunction, snm.obj, model.objects, M.matrix, lf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rff}{
 The rff slot from the lmer output
}
  \item{basisSplineFunction}{
B-spline basis function defined by buildBasisFunction
}
  \item{snm.obj}{
An object of class snm
}
  \item{model.objects}{
 A list containing the formatted model matrices returned from make.ref.model.matrices
}
  \item{M.matrix}{
 Matrix of pooled estimated RNA concentrations. The pooling strategy is defined by the nbins parameter. The element in position i,j is equal to the average RNA concentration on array j of the probes' whose average concentration is in the region spanned by bin i.
}
  \item{lf}{
 Object of class mer used for parsing random effects object.
}
}



\author{
Brig Mecham <brig.mecham@sagebase.org>
}




\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (rff, basisSplineFunction, snm.obj, model.objects, M.matrix) 
{
  splineDimNames <- paste("Bt", 1:snm.obj$spline.dim, sep = "")
  ranfx <- rff
  ranfx2 <- lapply(ranfx, function(x) {
    x[, splineDimNames]
  })
  ranfx2 <- ranfx2[names(snm.obj$int.var)]
  for (i in 1:length(snm.obj$int.var)) {
    x <- names(snm.obj$int.var)[i]
    ranfx2[[x]] <- as.matrix(ranfx2[[x]][levels(snm.obj$int.var[[x]]),])
  }
  ars <- sapply(1:dim(M.matrix)[2], function(i) {
    mREFs <- sapply(1:length(ranfx2), function(j) {
      model.objects$F.mats[[j]][i, ] \%*\% ranfx2[[j]]
    })
    bSM <- predict(basisSplineFunction, M.matrix[, as.numeric(i)])
    arsL <- bSM \%*\% mREFs
    rowSums(arsL)
  })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
