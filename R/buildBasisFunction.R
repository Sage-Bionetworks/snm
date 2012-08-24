
buildBasisFunction <-
		function(snm.obj) 
{
  lower <- min(rowMeans(snm.obj$dat))
  upper <- max(rowMeans(snm.obj$dat))
  basisSplineFunction <- ns(seq(lower,upper,length=100), df = snm.obj$spline.dim, 
      Boundary.knots = c(lower - 1, upper + 1),intercept=TRUE)
}
