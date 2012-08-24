make.snm.obj <-
		function(Y, bio.var, adj.var,int.var, spline.dim, nbins, weights=NULL, diagnose, rm.adj)
{
  if(!is.matrix(Y)) {
    stop("Y must be a matrix.")
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
    if(nrow(bio.var)!=ncol(Y)) {
      stop("The dimensions of Y and bio.var are incompatible. Please correct these objects.")
    }
    bio.var = cbind(bio.var[,!apply(bio.var, 2, function(x) {length(unique(x))==1})])
    if(ncol(bio.var)==0) {
      stop("bio.var must have at least one column (that is not an intercept term).")
    }
  }
  
  if(is.null(adj.var)) { 
    adj.var = cbind(rep(1,ncol(Y)))
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
  if(nrow(adj.var)!=ncol(Y)) {
    stop("The dimensions of Y and adj.var are incompatible. Please correct these objects.")
  }
  adj.var = cbind(rep(1,ncol(Y)), adj.var[,!apply(adj.var, 2, function(x) {length(unique(x))==1})])

  if(!is.null(int.var)) {
    if(is.vector(int.var)) {int.var = data.frame(as.factor(int.var))}
    if(is.matrix(int.var)) {
      int.var.temp = data.frame(as.factor(int.var[,1]))
      if(ncol(int.var)>1) {
        for(i in 2:ncol(int.var)) {
          int.var.temp = data.frame(int.var.temp, as.factor(int.var[,i]))
        }
      }
      if(!is.null(dimnames(int.var)[[2]])) {dimnames(int.var.temp)[[2]] = dimnames(int.var)[[2]]}
      int.var = int.var.temp
    }
    if(!is.data.frame(int.var)) {
      stop("int.var must be a data frame")
    }
    for(i in 1:ncol(int.var)) {
      if(!is.factor(int.var[,i])) {
        stop("Each column of int.var must be a factor variable.")
      }
    }
  }

  xx = try(solve(t(adj.var)%*%adj.var), silent=TRUE)
  if(is.character(xx)) {
    stop("adj.var is not a valid model matrix. Enter '?model.matrix' for more information on building a model matrix.\n")
  }
  yy = try(solve(t(cbind(bio.var,adj.var))%*%cbind(bio.var,adj.var)), silent=TRUE)
  if(is.character(yy)) {
    stop("cbind(bio.var,adj.var) is not a valid model matrix. Enter '?model.matrix' for more information on building a model matrix.\n")
  }
  
  dimnames.adj = dimnames(adj.var)
  dimnames(adj.var) = list(dimnames(adj.var)[[1]], paste("A", 1:ncol(adj.var), sep="")) 
  dimnames.bio = dimnames(bio.var)
  if(!is.null(bio.var)) {dimnames(bio.var) = list(dimnames(bio.var)[[1]], paste("B", 1:ncol(bio.var), sep=""))} 
  
  snm.obj = list()
  snm.obj$n.arrays = ncol(Y)
  snm.obj$n.probes = nrow(Y)
  snm.obj$rm.adj = rm.adj
  snm.obj$bio.var = bio.var
  snm.obj$adj.var = adj.var
  if(!is.null(bio.var)) {
    snm.obj$df.full = ncol(bio.var) + ncol(adj.var)
  } else {
    snm.obj$df.full = ncol(adj.var)
  }
  snm.obj$df.null = ncol(adj.var)
  snm.obj$int.var = int.var
  snm.obj$individuals = NULL
  snm.obj$spline.dim = spline.dim
  snm.obj$nbins = nbins
  snm.obj$diagnose = diagnose
  snm.obj$dat=Y
  if(is.null(weights)) {weights <- rep(1,nrow(Y))}
  snm.obj$weights=weights
  snm.obj$r.dat=Y
  snm.obj$dimnames.adj = dimnames.adj
  snm.obj$dimnames.bio = dimnames.bio
  class(snm.obj) = "edge"
  return(snm.obj)
}
