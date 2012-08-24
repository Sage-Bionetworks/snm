fit.pset <-
		function(dat, pset2row,downweight=TRUE, verbose=TRUE) { 
	# Sweep out the row means
	dat.c <- sweep(dat,1,rowMeans(dat))
	n <- ncol(dat)
	mns <- vars <- matrix(0, nr=length(pset2row), nc=n)
	singular.values <- matrix(0, nr=length(pset2row), nc=3)
	rownames(mns) <- rownames(vars) <- rownames(singular.values) <- names(pset2row)
	colnames(mns) <- colnames(vars) <- colnames(dat)
	colnames(singular.values) <- paste("Singular Value", 1:ncol(singular.values))
	probe.weights <- rep(0, nrow(dat))
	ptm <- proc.time()
	oneToOne <- TRUE
	if(sum(sapply(pset2row,length)) != length(probe.weights)){
		# This implies either a one-to-many mapping or the existing of some probes that don't map to a probe set.
		# Here we store the probe weights as a list.
		oneToOne <- FALSE
		probe.weights <- vector("list", length(pset2row))
		names(probe.weights) <- names(pset2row)
	}
# fit models for each probe
	for(i in 1:length(pset2row)) {		
		if((i %% 1000)==0 & verbose) {
			cat("\rCompleted ", i, " of ", length(pset2row), " total probe sets.")
		}
		rows <- pset2row[[i]]
		m <- length(rows)
		mat <- dat[rows,]
		if(length(rows) == 1){
			mns[i,] <- dat[rows,]
			vars[i,] <- rep(0, ncol(dat))
			if(!oneToOne){
				probe.weights[[names(pset2row)[i]]] <- 0
			}else{
				probe.weights[rows] <- 0
			}
			singular.values[i,] <- rep(0,3)
		}else if (length(rows) == 2){
			mns[i,] <- apply(mat,2,mean)
			vars[i,] <- apply(mat,2,var)
			if(!oneToOne){
				probe.weights[[names(pset2row)[i]]] <- c(summary(lm(mat[1,] ~ mns[i,]))$r.squared, 
						summary(lm(mat[2,] ~ mns[i,]))$r.squared)
			}else{
				probe.weights[rows] <- c(summary(lm(mat[1,] ~ mns[i,]))$r.squared, 
						summary(lm(mat[2,] ~ mns[i,]))$r.squared)
			}
			singular.values[i,] <- c(mean(probe.weights[rows]),0,0)
		}else{
			if(downweight) { 
				mat <- downweight.outliers(mat)
			}
			mat.c <- mat - rowMeans(mat)
			if(class(try(u <- fs(mat),silent=TRUE)) != "try-error") {
				if(sum(is.na(u$v[,1])) == 0) {
					X <- model.matrix(~u$v[,1])		
					cfs <- solve(t(X) %*% X) %*% t(X) %*% t(mat)
					est.data <- t(X %*% cfs)
					wts <- rowSums((mat - est.data)^2) / rowSums(mat.c^2)
					wts[wts>0.9999] <- 0.9999	
					
					cfs <- colSums(mat*(1-wts)) / sum(1-wts)
					res <-t(t(mat) - cfs)    
					res2 <- sweep(res,1,rowMeans(res))
					vr <- colSums((1-wts)*res2^2) / sum(1-wts)
					
					res.mat <- matrix(res2, nr=m)
					mns[i,] <- cfs
					vars[i,] <- vr
					if(!oneToOne){
						probe.weights[[names(pset2row)[i]]] <- wts
					}else{
						probe.weights[rows] <- wts
					}
					singular.values[i,] <- u$d[1:3]
				}else{
					mns[i,] <- apply(mat,2,mean)
					vars[i,] <- apply(mat,2,var)
					if(!oneToOne){
						probe.weights[[names(pset2row)[i]]] <- rep(0,length(rows))
					}else{
						probe.weights[rows] <- rep(0,length(rows))
					}
					singular.values[i,] <- rep(1,3)
				}
			}else{			
				mns[i,] <- apply(mat,2,mean)
				vars[i,] <- apply(mat,2,var)
				if(!oneToOne){
					probe.weights[[names(pset2row)[i]]] <- rep(0,length(rows))
				}else{
					probe.weights[rows] <- rep(0,length(rows))
				}
				singular.values[i,] <- rep(2,3)
			}
		}
	}
	if(oneToOne){
		names(probe.weights) <- rep(names(pset2row), sapply(pset2row,length))
		ret.obj <- list(estimated.rna.concentration = mns,
				estimated.rna.variances = vars,
				singular.values = singular.values,
				probe.weights = 1 - probe.weights)
	}else{
		ret.obj <- list(estimated.rna.concentration = mns,
				estimated.rna.variances = vars,
				singular.values = singular.values,
				probe.weights = lapply(probe.weights, function(x){ 1 - x }))
	}
	if(verbose){ cat("\n") }
	class(ret.obj) <- "snmPSET"
	ret.obj
}

