sim.probe.specific <- function(data=NULL,variable=NULL,rows=NULL,sample.from=NULL) {  
	
#  Make certain only a single variable was passed
  if(class(variable) != "character" & class(variable) != "numeric" & class(variable) != "factor") { 
    stop("Variable not valid.  Must be a character, vector or factor.")
  }
	
#Make any character a factor
  if(class(variable)=="character") {
    variable <- as.factor(variable)
  }
	
# Handle the rows variable.  If its a single number, then its a proportion.  Otherwise, its a vector of probe indices to be modified.
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
  
  x <- model.matrix(~-1+variable)  # Create model matrix for variable
  sample.this.many <- sum(lab==1) * ncol(x)  # Count how many probes to be modified
  cfs <- matrix(0, nr=length(lab), nc=ncol(x))  #Initialize matrix of coefficient
  # If sample.from is a list, sample from the function and pass the parameters.
  if(is.list(sample.from)) {
    #Estimate coefficients and add to coefficients matrix 
    cfs[which(lab==1),] <- do.call(sample.from$func, as.list(c(n=sample.this.many, sample.from$params)))
  }else{
    stop("sample.from must be a list")
  }
  cfs.mat <- cfs %*% t(x)  #  Estimate overall effect
  t(apply(cfs.mat,1,function(x) {  #Make all positive.  Don't want streaks in scatter plots
    				if(min(x) < 0) {
      				x - min(x)
    				}else{
      				x
    				}
  				})) -> cfs.mat
  cfs.mat  # Return effects
}
