get.coefs <-
		function(Y, mod) 
{
	t(solve(t(mod) %*% mod) %*% t(mod) %*% t(Y))
}
