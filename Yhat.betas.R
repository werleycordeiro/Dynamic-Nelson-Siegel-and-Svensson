Yhat.betas <- function(Y,Z)
{
  beta <- solve(t(Z) %*% Z) %*% t(Z) %*% t(Y) # Z=maturities X betas, Z'=betas X maturities, Y=observations X maturities, Y'=maturities X observations
  beta <- t(beta)
  Yhat <- beta %*% t(Z)
  
  return(list(Yhat=Yhat,beta=beta))
}
