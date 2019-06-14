NS <- function(Y)
{
  beta <- solve(t(Z) %*% Z) %*% t(Z) %*% t(Y) # Z=17x3, Z'=3x17, Y=348x17, Y'=17x348
  beta <- t(beta)
  Yhat <- beta %*% t(Z)
  
  return(list(Yhat=Yhat,beta=beta))
}
