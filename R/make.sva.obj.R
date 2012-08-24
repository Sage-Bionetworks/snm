make.sva.obj <-
function(dat, bio.var, adj.var, n.sv, num.iter, diagnose)
{		
  if(!is.matrix(dat)) {
    stop("dat must be a matrix.")
  }
	if(is.null(num.iter) | num.iter<1) { 
		stop("Please specify the number of iterations")
	}	
  if(!is.null(bio.var)) {
    if(is.data.frame(bio.var)) {
      bio.temp = NULL
      for(i in 1:ncol(bio.var)) {
        bio.temp = cbind(bio.temp, model.matrix(~bio.var[,i])[,-1])
      }
      bio.var = bio.temp
    }
    if(!is.matrix(bio.var)) {
      stop("bio.var must be a matrix.")
    }
    if(nrow(bio.var)!=ncol(dat)) {
      stop("The dimensions of Y and bio.var are incompatible. Please correct these objects.")
    }
    bio.var = cbind(bio.var[,!apply(bio.var, 2, function(x) {length(unique(x))==1})])
    if(ncol(bio.var)==0) {
      stop("bio.var must have at least one column (that is not an intercept term).")
    }
  }
  
  if(is.null(adj.var)) { 
    adj.var = cbind(rep(1,ncol(dat)))
  }
  if(is.data.frame(adj.var)) {
    adj.temp = NULL
    for(i in 1:ncol(adj.var)) {
      adj.temp = cbind(adj.temp, model.matrix(~adj.var[,i])[,-1])
    }
    adj.var = adj.temp
  }
  if(!is.matrix(adj.var)) {
    stop("adj.var must be a matrix.")
  }
  if(nrow(adj.var)!=ncol(dat)) {
    stop("The dimensions of Y and adj.var are incompatible. Please correct these objects.")
  }
  adj.var = cbind(rep(1,ncol(dat)), adj.var[,!apply(adj.var, 2, function(x) {length(unique(x))==1})])
	
  sva.obj = list()
	sva.obj$dat = dat
  sva.obj$n.arrays = ncol(dat)
  sva.obj$n.probes = nrow(dat)
  sva.obj$bio.var = bio.var
  sva.obj$adj.var = adj.var
	if(!is.null(bio.var)) {
  	sva.obj$df.full = ncol(bio.var) + ncol(adj.var)
	} else {
  	sva.obj$df.full = ncol(adj.var)
	}	
	sva.obj$df.null = ncol(adj.var)
	sva.obj$individuals = NULL
	sva.obj$num.iter=num.iter
	sva.obj$n.sv=n.sv
	sva.obj$svd = list()
	sva.obj$iteration = 1
	sva.obj$diagnose=diagnose
	sva.obj$singular.values = matrix(NA, nr=sva.obj$num.iter, nc=ncol(dat))
	class(sva.obj) = "edge"
  return(sva.obj)
}

