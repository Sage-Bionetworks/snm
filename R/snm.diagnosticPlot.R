snm.diagnostic.plot <-
  function(obs.fit, pval, snm.obj, iter)
{
  par(mfrow=c(2,2), oma=c(1,0,2,0))
  bds <- ifelse(snm.obj$bin.densities > 500, 500, snm.obj$bin.densities)
  plot(rowMeans(snm.obj$M.pooled),
       bds,ylim=c(0,500), pch=19, yaxt="n",
       xlab="Estimated RNA Concentration", ylab="Probes per bin", main="Null Probes Per Bin")
  axis(side=2, at=500, label=" 500+")
  axis(side=2, at=c(0,100,200,300,400))
  
  u <- fast.svd(obs.fit$res1,tol=0)
  plot(100*u$d^2 / sum(u$d^2), ylim=c(0,100), pch=19,
       main="Residual Latent Structure",
       xlab="Principal Component",
       ylab="Percent Variation Explained")
  
  if(dim(snm.obj$dat)[1] > 5000) {
    th <- sample(dim(snm.obj$dat)[1], 5000)
  }else{
    th <- 1:dim(snm.obj$dat)[1]
  } 
  oo <- order(snm.obj$M[th,1])
  plot(snm.obj$M[th[oo],1], snm.obj$array.fx[th[oo],1],
       type="l",lwd=1, ylim=range(snm.obj$array.fx[th[oo],]),
       xlab="Estimated Intensity", ylab="Estimated Effect by Array",main="Intensity-Dependent Effects")
  sapply(2:dim(snm.obj$r.dat)[2],function(id) {
    oo <- order(snm.obj$M[th,id])
    points(snm.obj$M[th[oo],id], snm.obj$array.fx[th[oo],id],type="l")
  }) -> hmm
  
  hist(pval, xlab="P-values",main="P-value Distribution");
  abline(v=min(pval[snm.obj$nulls]),col="red")
  pi0 = round(snm.obj$pi0,3)
  mtext(substitute(hat(pi)[0] == that, list(that=pi0)))
  title(paste("SNM Diagnostic Iteration", iter, sep=" "), outer=TRUE)
}
