calcArrayEffects <-
  	function(rff, basisSplineFunction, snm.obj, model.objects, M.matrix,lf) 
{
  offset <- 1;
  rfx <- list();
  for(i in 1:ncol(snm.obj$int.var)) { 
    rfx[[i]] <- matrix(rff@ranef[offset:(-1 + offset + nrow(lf$FL$trms[[i]]$Zt))], nr=length(unique(snm.obj$int.var[,i])))
    offset <- offset + nrow(lf$FL$trms[[i]]$Zt)
  }  
  ars <- sapply(1:ncol(M.matrix), function(i) {
    		mREFs <- sapply(1:length(rfx), function(j) {
      				model.objects$F.mats[[j]][i, ] %*% rfx[[j]]
    				})
    		bSM <- predict(basisSplineFunction, M.matrix[, as.numeric(i)])
    		arsL <- bSM %*% mREFs
    		rowSums(arsL)
  		})
  ars
}
