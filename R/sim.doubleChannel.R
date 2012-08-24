sim.doubleChannel <- function(seed) { 
  set.seed(seed)
  np <- 25000
  na <- 20
  gmeans <- rchisq(np,1,2)
  gmeans[gmeans>15] <- runif(sum(gmeans>15),15,16)
  data <- matrix(gmeans,nr=np,nc=na)
  bio.var <- data.frame(groups=rep(c("A","B"),each=10))
  adj.var <- data.frame(height=rnorm(20,1,0.5))
  int.var <- data.frame(array=factor(c(1:10,1:10)), dye=factor(rep(c("CY3","CY5","CY5","CY3"),each=5)))
  retBio<-FALSE
  
  group.effect <- sim.probe.specific(data, bio.var$groups, 0.3, list(func=rnorm,params=c(mean=1,sd=0.3)))
  height.effect <- sim.probe.specific(data, adj.var$height, 0.2, list(func=rnorm, params=c(mean=1,sd=0.1)))
  
  M <- data + group.effect + height.effect
  
  array.effect <- sim.intensity.dep(M, int.var$array, 2, list(func=rnorm, params=c(mean=0,sd=1)))
  dye.effect <- sim.intensity.dep(M, int.var$dye, 2, list(func=rnorm, params=c(mean=0,sd=1)))
  E <- matrix(rnorm(length(data),0,0.25), nr=nrow(data), nc=ncol(data))
  Y <- M + array.effect + dye.effect + E
  
  true.nulls <- which(group.effect[,1] == group.effect[,11])
  ret.obj <- 
    	list(raw.data=Y, 
          bio.var=model.matrix(~groups,data=bio.var),
          adj.var=model.matrix(~height,data=adj.var),
          int.var=int.var,
          true.nulls=true.nulls)
  ret.obj
}
