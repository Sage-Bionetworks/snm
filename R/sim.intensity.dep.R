sim.intensity.dep <- function(data=NULL,variable=NULL,spline.dim=NULL,sample.from=NULL, rows=1) {
	
#  Make certain only a single variable was passed
  if(class(variable) != "character" & class(variable) != "factor") { 
    stop("Variable not valid.  Must be a character or factor.")
  }
  
#Make any character a factor
  if(class(variable)=="character") {
    variable <- as.factor(variable)
  }
  
# Handle the rows variable.  If its a single number, then its a proportion.  Otherwise, its a vector of probe indices to be modified.
  if(rows!=1) { 
    lab <- rep(0, nrow(data))  # Create vector of dummy variables    
    if(length(rows)==1) {
      if(rows > 1) { stop("Proportion of data to be influenced by variable greater than 1.")}
      these <- sample(1:nrow(data), nrow(data) * rows) # Sample rows to be modified
      lab[these] <- 1  #Identify rows to be modified
    }else{
      if(max(rows) > nrow(data)) { stop("At least one row of data matrix to be influenced by variable is larger than total number of rows in data")}
      if(min(rows) < 0) { stop("Rows to be influenced by variable must be greater than 0")}
      lab[rows] <- 1  # Identify rows to be modified
    }
  }else{
    lab=rep(1,nrow(data))
  }
  
  x <- model.matrix(~-1+variable)  # Create model matrix for variable
  z <- ns(seq(min(data),max(data),length=100), df=spline.dim, Boundary.knots=c(min(data)-1,max(data)+1))
  cfs.mat <- 0 * data
  spline.coefficients <- matrix(0,nr=ncol(x), nc=spline.dim)
  for(i in 1:ncol(x)) { 
    spline.coefficients[i,] <- do.call(sample.from$func, as.list(c(n=ncol(z), sample.from$params)))
  }
  for(i in 1:ncol(data)) { 
    spline.evals <- predict(z,data[,i])
    cfs.mat[,i] <- spline.evals %*% spline.coefficients[which(x[i,]==1),]
  }
  cfs.mat <- lab * cfs.mat
  cfs.mat  # Return coefficients
}
