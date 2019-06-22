VARcoeff<-function(betas,l)
{
  YY <-betas
  XX <- cbind(rep(1,length(betas[,1])),YY)
  XX <-XX[-l,] # exclude last observation
  YY <-YY[-1,] # exclude first observation
  var <- solve(t(XX) %*% XX) %*% t(XX) %*% YY
  var <-t(var)
  
  var # Start point to pars$phi in Kalman-Filter-Dynamic-Nelson-Siegel
}
