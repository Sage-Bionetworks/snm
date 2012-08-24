
calcSig <-
		function(Y=NULL, bio.var=NULL, adj.var=NULL, df.full=NULL, df.null=NULL) { 
	if(!is.matrix(Y)) {
		stop("Y must be a matrix.")
	}
	if(is.null(adj.var)) { 
		adj.var = cbind(rep(1,ncol(Y)))
	} else if(!is.matrix(adj.var)) {
		stop("adj.var must be a matrix.")
	}
	if(is.null(bio.var)) {
		stop("bio.var must be supplied by the user (and cannot be an intercept term only).")
	} else if(!is.matrix(bio.var)) {
		stop("bio.var must be a matrix.")
	}
	if(!all(nrow(bio.var)==nrow(adj.var), nrow(bio.var)==ncol(Y))) {
		stop("The dimensions of Y, bio.var, and adj.var are incompatible. Please correct these objects.")
	}
	if(!all(is.matrix(bio.var), is.matrix(adj.var), is.matrix(Y))) {
		stop("bio.var, adj.var, and Y must all be matrices.")
	}
	if(!is.null(bio.var)) {
		bio.var = cbind(bio.var[,!apply(bio.var, 2, function(x) {length(unique(x))==1})])
		if(ncol(bio.var)==0) {
			stop("bio.var must have at least one column (that is not an intercept term).")
		}
	}
	adj.var = cbind(rep(1,ncol(Y)), adj.var[,!apply(adj.var, 2, function(x) {length(unique(x))==1})])
	inf.obj = list()
	inf.obj$n.arrays = ncol(Y)
	inf.obj$n.probes = nrow(Y)
	inf.obj$bio.var = bio.var
	inf.obj$adj.var = adj.var
	if(is.null(df.full)) {
		inf.obj$df.full = ncol(bio.var) + ncol(adj.var)
	}else{
		inf.obj$df.full = df.full
	}
	if(is.null(df.null)) { 
		inf.obj$df.null = ncol(adj.var)
	}else{
		inf.obj$df.null = df.null
	}
	inf.obj$dat = Y
	class(inf.obj) = "edge"
	
	obs.fit <- edge.fit(inf.obj, odp=FALSE)
	obs.stat <- edge.glr(obs.fit, df1=inf.obj$df.full, df0=inf.obj$df.null, norm.pval=TRUE)
	inf.obj$cfs <- t(obs.fit$c1)
	inf.obj$pval <- obs.stat$pval; 
	tmp <- try(edge.qvalue(obs.stat$pval))
	if(class(tmp) == "try-error"){
		inf.obj$pi0 <- NA
		inf.obj$qvalues <- NA
	}else{
		inf.obj$pi0 <- tmp$pi0
		inf.obj$qvalues <- tmp$qvalues
	}
	inf.obj$total.ssq = sum(sweep(Y,1,rowMeans(Y))^2)
	inf.obj$full.resid.ssq <- obs.stat$s.rss1
	inf.obj$red.resid.ssq <- obs.stat$s.rss0
	inf.obj$stat <- obs.stat$stat
	inf.obj$dat = NULL
	
	inf.obj
}

