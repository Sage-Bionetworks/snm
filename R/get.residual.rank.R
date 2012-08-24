get.residual.rank <-
function(sva.obj,...) {
  mod <- cbind(sva.obj$adj.var, sva.obj$bio.var)
	
	dat <- sva.obj$dat
  H <- mod %*% solve(t(mod) %*% mod) %*% t(mod)
  res <- dat - t(H %*% t(dat))
	B = crossprod(res)
	s = svd(B, nu = 0,nv=0)
	ds <- s$d / sum(s$d)
	if(is.null(sva.obj$n.sv)) {
	p.values <- sapply(1:(ncol(dat)-ncol(mod)-5), function(x) {
				ks.test(pnorm(ds[x:(ncol(dat)-ncol(mod))], 
								mean(ds[x:(ncol(dat)-ncol(mod))]), 
								sd(ds[x:(ncol(dat)-ncol(mod))])),"punif")$p})
	estimated.rank <- sum(p.values < 0.05)
	sva.obj$n.sv = estimated.rank
	}else{
		estimated.rank = sva.obj$n.sv
	}
	s = svd(B, nu = 0, nv=estimated.rank)
	s$d <- s$d / sum(s$d)	

	sva.obj$svd[[1]] = s
	sva.obj$singular.values[1,] <- s$d
  cat("Estimated rank of dependence kernel: ", estimated.rank,"\n")
	sva.obj
}

