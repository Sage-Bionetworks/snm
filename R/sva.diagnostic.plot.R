sva.diagnostic.plot <-
function(ptmp.b, pprob.b, ptmp.u, pprob.u, pprob, sva.obj)
{
  par(mfrow=c(2,2), oma=c(1,0,2,0))	
	tmp <- hist(ptmp.b,plot=FALSE)
	tmp$counts <- tmp$counts/sum(tmp$counts)
	plot(tmp,xlab="P-values",main="Observed Influence")
	oo <- order(ptmp.b)
	points(ptmp.b[oo], pprob.b[oo],col="red",type="l",lwd=3)
	
	tmp <- hist(ptmp.u,plot=FALSE)
	tmp$counts <- tmp$counts/sum(tmp$counts)
	plot(tmp,xlab="P-values",main="Latent Influence")
	oo <- order(ptmp.u)
	points(ptmp.u[oo], pprob.u[oo],col="red",type="l",lwd=3)
	
	tmp <- hist(pprob,plot=FALSE)
	tmp$counts <- tmp$counts/sum(tmp$counts)
	plot(tmp,xlab="Probabilities",main="Probabilities used \nto weight SVD")
	
	iter <- length(sva.obj$svd)
	plot(sva.obj$singular.values[,1], 
			type="l", 
			main="Singular Values\nby iterations", 
			xlab="Iteration", 
			ylab="Singular Values or Eigenweights", 
			ylim=c(0, max(sva.obj$singular.values[,1],na.rm=TRUE)))
	for(i in 2:sva.obj$n.sv) {
		points(sva.obj$singular.values[,i], type="l")
	}
}