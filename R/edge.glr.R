edge.glr <-
function(obs.fit, df1, df0, norm.pval=FALSE) {
  #Function by John Storey for the edge package
  err.func = "edge.glr"
	
  res1 = obs.fit$res1
  res0 = obs.fit$res0
  
  n = ncol(res1)
  
  rss1 = apply(res1, 1, function(x) {sum(x^2)})
  rss0 = apply(res0, 1, function(x) {sum(x^2)})
	
  stat = ((rss0 - rss1)/(df1-df0))/(rss1/(n-df1))
	
  if(norm.pval) {
    pval = 1 - pf(stat, df1=(df1-df0), df2=(n-df1))
    return(list(stat=stat, pval=pval, s.rss1=sum(rss1), s.rss0=sum(rss0)))
  } else {
    return(list(stat=stat, pval=NULL))
  }
	
}

