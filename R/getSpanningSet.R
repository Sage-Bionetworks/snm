getSpanningSet <-
		function(snm.obj) {
  m <- rowMeans(snm.obj$M[snm.obj$nulls,])
  steps <- seq(quantile(snm.obj$M[snm.obj$nulls,],0.001),
      quantile(snm.obj$M[snm.obj$nulls,],0.999), length.out=snm.obj$nbins)
  th <- sapply(m, function(x) { 
    		sum(x <= steps)
  		}) 
  bins <- split(snm.obj$nulls, th)
  bins
}