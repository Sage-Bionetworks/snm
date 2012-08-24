sim.refDesign <- function(seed) { 
  set.seed(seed)
  np <- 25000
	na <- 40
	gmeans <- rchisq(np,1,2)
	gmeans[gmeans>15] <- runif(sum(gmeans>15),15,16)
	data <- matrix(gmeans,nr=np,nc=na)
	bio.var <- data.frame(groups=rep(c("A","B","C"),c(10,10,20)))
	adj.var <- NULL
	int.var <- data.frame(array=factor(c(1:20,1:20)), dye=factor(rep(c("CY3","CY5"),each=20)))
	
	group.effect <- sim.probe.specific(data, bio.var$groups, 0.3, list(func=rnorm,params=c(mean=1,sd=0.3)))
	
	M <- data + group.effect 
	
	array.effect <- sim.intensity.dep(M, int.var$array, 2, list(func=rnorm, params=c(mean=0,sd=1)))
	dye.effect <- sim.intensity.dep(M, int.var$dye, 2, list(func=rnorm, params=c(mean=0,sd=1)))
	E <- matrix(rnorm(length(M),0,0.25), nr=nrow(M), nc=ncol(M))
	Y <- M + array.effect + dye.effect + E
	
	true.nulls <- which(group.effect[,1] == group.effect[,11])
	ret.obj <- 
			list(raw.data=Y, 
					bio.var=model.matrix(~groups,data=bio.var),
					adj.var=NULL,
					int.var=int.var,
					true.nulls=true.nulls)
	ret.obj
}
