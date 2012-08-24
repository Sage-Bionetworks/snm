buildBasisSplineMatrix <-
  	function(M, basisSplineFunction) 
{
  lnp <- dim(M)[1]
  bSM <- apply(M, 2, function(x) {
    		predict(basisSplineFunction, x)
  		})
  bSM.model <- sapply(1:(dim(basisSplineFunction)[2]), function(x) {
    		pos <- (x - 1) * lnp + 1
    		as.numeric(bSM[pos:(pos + lnp - 1), ])
  		})
  colnames(bSM.model) <- paste("Bt", 1:dim(bSM.model)[2], sep = "")
  bSM.model
}