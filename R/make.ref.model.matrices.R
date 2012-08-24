make.ref.model.matrices <-
		function(snm.obj, exp) 
{
  F.mats <- list()
  spline.re.model <- ""
  for (i in 1:length(snm.obj$int.var)) {
    ide <- names(snm.obj$int.var)[i]
    spline.re.model <- paste(spline.re.model, "+ ( sp -1 | ", ide, ")", sep = "")
    F.mats[[i]] <- model.matrix(~-1 + snm.obj$int.var[,i])
  }
  if (dim(snm.obj$adj.var)[2] > 0) {
    ref.model <- paste( "probes+", paste("probes:", colnames(snm.obj$adj.var)[-1],collapse = "+"),sep="+" )
  }
  else {
    ref.model <- "probes + "
  }
  ref.model <- as.formula(paste("y ~ -1 +", ref.model, spline.re.model), 
      env = exp)
  names(F.mats) <- names(snm.obj$int.var)
  list(ZF = ref.model, F.mats = F.mats)
}
