\name{fit.pset}
\alias{fit.pset}
\title{
Summarize probe-level data into probe sets
}
\description{
This function builds probe-set level summaries for probe-level data. Requires a matrix of probe-level data,  
and a list describing the probe to probe set mappings.
}
\usage{
fit.pset(dat, pset2row, downweight = TRUE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
An \eqn{m} probes by \eqn{n} arrays matrix of expression data. 
}
  \item{pset2row}{
A list with element \eqn{i} containg the indices in \eqn{dat} that correspond to the \eqn{i}th probe set.
}
  \item{downweight}{
A logical scalar.  If 'TRUE' outlier measurements are downweighted by the algorithm.
}
  \item{verbose}{
A logical scalar.  If ÔTRUEÕ verbose output is generated during the probe set summarization step.
}
}
\details{
This function implements the probe set summarization algorithm utilized by the snm framework.
}
\value{
  \item{estimated.rna.concentration }{
The matrix of summarized probe set values.
}
  \item{estimated.rna.variances}{
The matrix of probe set variances associated with each probe set and sample.
}
  \item{singular.values}{
An \eqn{g} probe sets by \eqn{3} matrix consisting of the first three singular values obtained from the decomposition
of each probe sets probe-level data. 
}
  \item{probe.weights}{
A vector of length \eqn{m} containing the probe weights.  Each weight describes the relationship between the data for a given probe and its summarized probe set vector.  The weights range from [0,1], with weights close to 0 implying stronger relationships.
}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brig Mecham <brig.mecham@sagebase.org> 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dat, pset2row, downweight = TRUE, verbose = TRUE) 
{
    dat.c <- sweep(dat, 1, rowMeans(dat))
    n <- ncol(dat)
    mns <- vars <- matrix(0, nr = length(pset2row), nc = n)
    singular.values <- matrix(0, nr = length(pset2row), nc = 3)
    rownames(mns) <- rownames(vars) <- rownames(singular.values) <- names(pset2row)
    colnames(mns) <- colnames(vars) <- colnames(dat)
    colnames(singular.values) <- paste("Singular Value", 1:ncol(singular.values))
    probe.weights <- rep(0, nrow(dat))
    ptm <- proc.time()
    for (i in 1:length(pset2row)) {
        if ((i\%\%1000) == 0 & verbose) {
            cat("\rCompleted ", i, " of ", length(pset2row), 
                " total probe sets.")
        }
        rows <- pset2row[[i]]
        m <- length(rows)
        mat <- dat[rows, ]
        if (downweight) {
            mat <- downweight.outliers(mat)
        }
        mat.c <- mat - rowMeans(mat)
        if (class(try(u <- fs(mat), silent = TRUE)) != "try-error") {
            if (sum(is.na(u$v[, 1])) == 0) {
                X <- model.matrix(~u$v[, 1])
                cfs <- solve(t(X) \%*\% X) \%*\% t(X) \%*\% t(mat)
                est.data <- t(X \%*\% cfs)
                wts <- rowSums((mat - est.data)^2)/rowSums(mat.c^2)
                wts[wts > 0.9999] <- 0.9999
                lmerFit = TRUE
                lf <- ""
                cfs <- colSums(mat * (1 - wts))/sum(1 - wts)
                res <- t(t(mat) - cfs)
                res2 <- sweep(res, 1, rowMeans(res))
                vr <- colSums((1 - wts) * res2^2)/sum(1 - wts)
                res.mat <- matrix(res2, nr = m)
                mns[i, ] <- cfs
                vars[i, ] <- vr
                probe.weights[rows] <- wts
                singular.values[i, ] <- u$d[1:3]
            }
            else {
                mns[i, ] <- apply(mat, 2, mean)
                vars[i, ] <- apply(mat, 2, var)
                probe.weights[rows] <- rep(0, length(rows))
                singular.values[i, ] <- rep(1, 3)
            }
        }
        else {
            mns[i, ] <- apply(mat, 2, mean)
            vars[i, ] <- apply(mat, 2, var)
            probe.weights[rows] <- rep(0, length(rows))
            singular.values[i, ] <- rep(2, 3)
        }
    }
    ret.obj <- list(estimated.rna.concentration = mns, estimated.rna.variances = vars, 
        singular.values = singular.values, probe.weights = probe.weights, 
        meta.data = meta.data, platform = platform)
    if (!is.null(varName)) {
        a <- paste("fits", varName, sep = ".")
        assign(a, ret.obj)
        save(list = paste("fits", varName, sep = "."), file = paste("fits", 
            varName, "Rda", sep = "."))
    }
    if (verbose) {
        cat("\n")
    }
    class(ret.obj) <- "snmPSET"
    ret.obj
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
