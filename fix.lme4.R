# TODO: Add comment
# 
# Author: brigmecham
###############################################################################

set.seed(12346)
np <- 25000
na <- 50
gmeans <- rchisq(np,1,2)
gmeans[gmeans>15] <- runif(sum(gmeans>15),15,16)
data <- matrix(gmeans,nr=np,nc=na)
bio.var <- data.frame(groups=rep(c("A","B"),each=25))
adj.var <- data.frame(batches=rep(c("A","B","C","D","E"),times=10),
		height=rnorm(50,1,0.5))
int.var <- data.frame(array=factor(1:50))
retBio<-FALSE

group.effect <- sim.probe.specific(data, bio.var$groups, 0.3, list(func=rnorm,params=c(mean=1,sd=0.3)))
batches.effect <- sim.probe.specific(data, adj.var$batches, 0.1, list(func=rnorm,params=c(mean=0,sd=0.3)))
height.effect <- sim.function.var(data, adj.var$height, 0.2, list(list(func=rnorm, params=c(mean=1,sd=0.1)), list(func=rnorm, params=c(mean=0.3,sd=0.1))), 2) 

M <- data + group.effect + batches.effect + height.effect

array.effect <- sim.intensity.dep(M, int.var$array, 2, list(func=rnorm, params=c(mean=0,sd=1)))
E <- matrix(rnorm(length(data),0,0.25), nr=nrow(data), nc=ncol(data))
Y <- M + array.effect + E

true.nulls <- which(group.effect[,1] == group.effect[,26])

bio.var <- model.matrix(~groups, data=bio.var)
adj.var <- model.matrix(~batches+height,data=adj.var)
weights <- NULL
spline.dim=4
num.iter=20
nbins=50
verbose=TRUE
diagnose=TRUE
alg.type="em"

if(is.null(weights)) { weights <- rep(1,dim(Y)[1])}
snm.obj <- make.snm.obj(Y, bio.var, adj.var, int.var, spline.dim,
		nbins,weights,diagnose,alg.type,retBio)  
basisSplineFunction <- buildBasisFunction(snm.obj)
pi0s <- rep(0,num.iter)
obs.fit <- edge.fit(snm.obj, odp=FALSE)
obs.stat <- edge.glr(obs.fit, df1=snm.obj$df.full, df0=snm.obj$df.null, norm.pval=TRUE)
snm.obj$pi0 <- edge.qvalue(obs.stat$pval)$pi0
snm.obj$nulls <- true.nulls


snm.obj$M <- snm.obj$dat - obs.fit$res1
snm.obj$M[snm.obj$nulls,] <- obs.fit$fit0[snm.obj$nulls,]

# Split the data into nbins bins based on their mean intensities
bins <- getSpanningSet(snm.obj)

# Build the matrix of weighted raw data and matrix of weighted fitted values for each bin.
lnp <- length(bins)
np <- 1:lnp
Y.pooled <- 0*snm.obj$dat[np,]
M.pooled <- 0*snm.obj$M[np,]
for(i in 1:lnp) {
	Y.pooled[i,] = apply(matrix(snm.obj$r.dat[as.vector(bins[[i]]),], ncol=ncol(snm.obj$dat)),2,
			weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
	M.pooled[i,] = apply(matrix(snm.obj$M[as.vector(bins[[i]]),], ncol=ncol(snm.obj$M)),2,
			weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
}

# Build the basis spline matrix for the pooled coefficients.
bSM.model <- buildBasisSplineMatrix(M.pooled, basisSplineFunction)
exp <- new.env()
# Build the data object and fit the mixed effects model
expObj <- makeDataObject(Y.pooled, np, snm.obj,
		bSM.model, exp,bins)

batches <- rep("A",dim(expObj)[1])
batches[which(expObj$batchesB==1)] <- "B"
batches[which(expObj$batchesC==1)] <- "C"
batches[which(expObj$batchesD==1)] <- "D"
batches[which(expObj$batchesE==1)] <- "E"
batches <- as.factor(batches)

expObj$batches <- batches
expObj$sp <- as.matrix(expObj[,9:12])

lf <- lmer(y ~ probes+probes:batches + probes:height +
				(-1+sp | array),
		data=expObj,weights=weights,doFit=FALSE)
lf$FL$trms[[1]]$ST <- matrix(0,nr=1,nc=1)
rownames(lf$FL$trms[[1]]$ST) <- colnames(lf$FL$trms[[1]]$ST) <- "spline" 

rff <- do.call(lme4:::lmer_finalize,lf)

# Add useful variables to snm.obj
snm.obj$E.pooled <- matrix(rff@resid, nr=dim(Y.pooled)[1])
snm.obj$Y.pooled <- Y.pooled
snm.obj$M.pooled <- M.pooled
snm.obj$rff <- rff@ranef
snm.obj$bin.densities <- sapply(bins,length)

ranfx <- matrix(rff@ranef,nr=50)
M.matrix <- snm.obj$M
ars <- sapply(1:dim(M.matrix)[2], function(i) {
			bSM <- predict(basisSplineFunction, M.matrix[, as.numeric(i)])
			arsL <- bSM %*% ranfx[i,] 
			rowSums(arsL)
		})

snm.obj$dat <- snm.obj$r.dat - ars
obs.fit <- edge.fit(snm.obj, odp=FALSE)
obs.stat <- edge.glr(obs.fit, df1=snm.obj$df.full, df0=snm.obj$df.null, norm.pval=TRUE)
snm.obj$pi0 <- edge.qvalue(obs.stat$pval)$pi0
pval <- ks.test(obs.stat$pval[true.nulls],"punif")$p

write.table(pval,file="pvalues.R",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE)


#nt - Number of random effects terms in the model      
#n  - Number of observations
#p  - Number of fixed effects
#q  - Number of random effects
#s  - Number of parameters in nonlinear model function

lf <- lmer(y ~ probes+probes:batches + probes:height +
				(-1+Bt1 | array) + (-1+Bt2 | array) + (-1+Bt3 | array) + (-1+Bt4 | array),
		data=expObj,weights=weights,doFit=FALSE)


Z <- Matrix(t(model.matrix(~ -1 + Bt1:array + Bt2:array + Bt3:array + Bt4:array, data=expObj)))
A <- rBind(lf$FL$trms[[1]]$A, lf$FL$trms[[2]]$A, lf$FL$trms[[3]]$A, lf$FL$trms[[4]]$A)
lf$FL$trms[[1]]$Zt <- Z
lf$FL$trms[[1]]$A <- A

lf$FL$trms <- list(lf$FL$trms[[1]])
lf$FL$dims[1] <- 1
class(lf$FL$dims) <- "integer"

attr(lf$FL$fl,"assign") <- 1; 
mode(attr(lf$FL$fl,"assign")) <- "integer"

FL <- lf$FL

rff <- do.call(lme4:::lmer_finalize,lf)




set.seed(12346)
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

bio.var <- model.matrix(~groups, data=bio.var)
adj.var <- model.matrix(~height,data=adj.var)
weights <- NULL
spline.dim=4
num.iter=20
nbins=30
verbose=TRUE
diagnose=TRUE
alg.type="em"

if(is.null(weights)) { weights <- rep(1,dim(Y)[1])}
snm.obj <- make.snm.obj(Y, bio.var, adj.var, int.var, spline.dim,
		nbins,weights,diagnose,alg.type,retBio)  
basisSplineFunction <- buildBasisFunction(snm.obj)
pi0s <- rep(0,num.iter)
obs.fit <- edge.fit(snm.obj, odp=FALSE)
obs.stat <- edge.glr(obs.fit, df1=snm.obj$df.full, df0=snm.obj$df.null, norm.pval=TRUE)
snm.obj$pi0 <- edge.qvalue(obs.stat$pval)$pi0
snm.obj$nulls <- true.nulls

snm.obj$M <- snm.obj$dat - obs.fit$res1
snm.obj$M[snm.obj$nulls,] <- obs.fit$fit0[snm.obj$nulls,]

# Split the data into nbins bins based on their mean intensities
bins <- getSpanningSet(snm.obj)

# Build the matrix of weighted raw data and matrix of weighted fitted values for each bin.
lnp <- length(bins)
np <- 1:lnp
Y.pooled <- 0*snm.obj$dat[np,]
M.pooled <- 0*snm.obj$M[np,]
for(i in 1:lnp) {
	Y.pooled[i,] = apply(matrix(snm.obj$r.dat[as.vector(bins[[i]]),], ncol=ncol(snm.obj$dat)),2,
			weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
	M.pooled[i,] = apply(matrix(snm.obj$M[as.vector(bins[[i]]),], ncol=ncol(snm.obj$M)),2,
			weighted.mean, w=snm.obj$weights[as.vector(bins[[i]])])
}

# Build the basis spline matrix for the pooled coefficients.
bSM.model <- buildBasisSplineMatrix(M.pooled, basisSplineFunction)
exp <- new.env()
# Build the data object and fit the mixed effects model
expObj <- makeDataObject(Y.pooled, np, snm.obj, exp,bins)
expObj$sp <- as.matrix(bSM.model)
model.objects <- make.ref.model.matrices(snm.obj, exp)

lf <- do.call("lmer", list(model.objects$ZF, expObj, NULL,TRUE,list(),NULL,FALSE, FALSE,1:nrow(expObj),expObj$weights))

for(i in 1:ncol(int.var)) { 
	lf$FL$trms[[i]]$ST <- matrix(0,nr=1,nc=1)
	rownames(lf$FL$trms[[i]]$ST) <- colnames(lf$FL$trms[[i]]$ST) <- paste("spline",i,sep="")
}

rff <- do.call(lme4:::lmer_finalize,lf)

# Add useful variables to snm.obj
snm.obj$E.pooled <- matrix(rff@resid, nr=dim(Y.pooled)[1])
snm.obj$Y.pooled <- Y.pooled
snm.obj$M.pooled <- M.pooled
snm.obj$rff <- rff@ranef
snm.obj$bin.densities <- sapply(bins,length)

snm.obj$array.fx <- calcArrayEffects(rff, basisSplineFunction, snm.obj,model.objects,snm.obj$M)

snm.obj$dat <- snm.obj$r.dat - snm.obj$array.fx
obs.fit <- edge.fit(snm.obj, odp=FALSE)
obs.stat <- edge.glr(obs.fit, df1=snm.obj$df.full, df0=snm.obj$df.null, norm.pval=TRUE)
snm.obj$pi0 <- edge.qvalue(obs.stat$pval)$pi0
pval <- ks.test(obs.stat$pval[true.nulls],"punif")$p


