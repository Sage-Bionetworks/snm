checkArguments <-
  	function(edge.obj) 
{
  m <- dim(edge.obj$dat)[1]
  n <- dim(edge.obj$dat)[2]
  d <- dim(edge.obj$bio.var)[2]
  n.bc <- dim(edge.obj$bio.var)[1]
#    if (n.bc != n) {
#        print("Number of columns of Y and number of rows of bio.vars are not equal\n ")
#        stop()
#    }
#    n.r.c <- dim(edge.obj$adj.var)[1]
#    if (n.r.c != n && dim(edge.obj$adj.var)[1] != 0) {
#        print("Number of columns of Y and number of rows of technical.covariates are not equal\n ")
#        stop()
#    }
#    if (!is.matrix(edge.obj$dat)) {
#        print("Y is not a numeric matrix\n")
#        stop()
#    }
  return(edge.obj)    
}