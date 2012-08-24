
seeds <- floor(runif(100, 10,10000))

sapply(seeds, function(seed){ 
set.seed(seed)
np <- 25000
na <- 60
gmeans <- rchisq(np,1,2)
gmeans[gmeans>15] <- runif(sum(gmeans>15),15,16)
data <- matrix(gmeans,nr=np,nc=na)
bio.var <- data.frame(groups=rep(c("A","B"),each=30))
adj.var <- data.frame(batches=rep(c("A","B","C"),times=20),
    									gene.fusion=rep(c("A","B"),c(45,15)))
int.var <- data.frame(array=factor(1:na))

group.effect <- sim.probe.specific(data, bio.var$groups, 0.3, list(func=rnorm,params=c(mean=1,sd=0.3)))
batches.effect <- sim.probe.specific(data, adj.var$batches, 0.1, list(func=rnorm,params=c(mean=0,sd=0.3)))
gene.fusion.effect <- sim.probe.specific(data, adj.var$gene.fusion, 0.2, list(func=rnorm, params=c(mean=1,sd=0.1)))

M <- data + group.effect + batches.effect + gene.fusion.effect

array.effect <- sim.intensity.dep(M, int.var$array, 2, list(func=rnorm, params=c(mean=0,sd=1)))
E <- matrix(rnorm(length(data),0,0.25), nr=nrow(data), nc=ncol(data))
Y <- M + array.effect + E
true.nulls <- which(group.effect[,1] == group.effect[,31])

data.obj <- 
		list(raw.data=Y, 
				bio.var=model.matrix(~groups,data=bio.var),
				adj.var=model.matrix(~batches+gene.fusion,data=adj.var),
				int.var=int.var,
				true.nulls=true.nulls)

#
#  Fit full model to data
#
		
snm.fit.result.pure <- snm(data.obj$raw.data,
											data.obj$bio.var,
											data.obj$adj.var,
											data.obj$int.var,diagnose=FALSE)
ks.test(snm.fit.result.pure$pval[data.obj$true.nulls],"punif")$p
gc()
}) -> res

#
#  Fit model lacking adjustment variables
#
snm.fit.result <- snm(data.obj$raw.data,
											data.obj$bio.var,
											NULL,
											data.obj$int.var,num.iter=4)
ks.test(snm.fit.result$pval[data.obj$true.nulls],"punif")

#
# Implement SVA on this fit
#
library(sva)
svaobj <- sva(snm.fit.result$dat, data.obj$bio.var, cbind(rep(1,60)))

#
#Use Pr(u_i=0) as weights 
#
snm.fit.result2 <- snm(data.obj$raw.data,
		data.obj$bio.var,
		NULL,
		data.obj$int.var,num.iter=10, weights=1-svaobj$pprob.gam)
ks.test(snm.fit.result2$pval[data.obj$true.nulls],"punif")

#
#Use Estimated SVA on data 
#
uu <- svaobj$sv
obj <- data.frame(b1=uu[,1], b2=uu[,2], b3=uu[,3])
adj.var <- model.matrix(~ b1+b2+b3, data=obj)
snm.fit.result2b <- snm(data.obj$raw.data,
		data.obj$bio.var,
		adj.var,
		data.obj$int.var,num.iter=10, weights=1-svaobj$pprob.gam)
ks.test(snm.fit.result2b$pval[data.obj$true.nulls],"punif")


#
# Implement SVA on snm.fit.result2
#
svaobj2 <- sva(snm.fit.result2$dat, data.obj$bio.var, cbind(rep(1,60)))

#
#Fit model with updated dependence kernel
#
uu <- svaobj2$sv
obj <- data.frame(b1=uu[,1], b2=uu[,2], b3=uu[,3])
adj.var <- model.matrix(~ b1+b2+b3, data=obj)
snm.fit.result3 <- snm(data.obj$raw.data,
									 		data.obj$bio.var,
											adj.var,
											data.obj$int.var,num.iter=10)
ks.test(snm.fit.result3$pval[data.obj$true.nulls],"punif")

hist(snm.fit.result3$pval[data.obj$true.nulls])


