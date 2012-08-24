f.pvalue <-
function(dat,mod,mod0){
  # This is a function for performing
  # parametric f-tests on the data matrix
  # dat comparing the null model mod0
  # to the alternative model mod. 
  n <- ncol(dat)
  m <- nrow(dat)
  df1 <- ncol(mod)
  df0 <- ncol(mod0)
  p <- rep(0,m)
  Id <- diag(n)
  
  res <- dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*% t(mod))
  res0 <- dat %*% (Id - mod0 %*% solve(t(mod0) %*% mod0) %*% t(mod0))
  
  rss1 <- res^2 %*% rep(1,n)
  rss0 <- res0^2 %*% rep(1,n)
  
  fstats <- ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
  p <-  1-pf(fstats,df1=(df1-df0),df2=(n-df1))
  return(p)
}
