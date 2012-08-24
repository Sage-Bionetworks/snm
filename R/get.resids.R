get.resids <-
function (Y, mod) 
{
  Y - t(mod %*% solve(t(mod) %*% mod) %*% t(mod) %*% t(Y))
}

