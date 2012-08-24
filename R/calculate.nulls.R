calculate.nulls <-
  	function(pvals, pi0) {
  which(rank(pvals) > (1-pi0)*length(pvals))
}
