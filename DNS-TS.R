# author: Werley Cordeiro
# werleycordeiro@gmail.com
# test from Rstudio 2
# Packages
list.of.packages <- c("highfrequency","vars","YieldCurve")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(YieldCurve)
library(vars)
library(highfrequency)

# Data
data <-read.csv("https://www.dropbox.com/s/inpnlugzkddp42q/bonds.csv?dl=1",header = TRUE, sep = ";")
require(xts)
data <- as.matrix(data[,which(names(data)=="M3"):which(names(data)=="M120")])
datas <- seq(as.Date("1972/1/1"), by = "month", length.out = 348)
data  <- xts(data, order.by = datas)
head(data)
dim(data)

l <-dim(data)[1] # observations

maturity<-c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120) # 17 maturities in the data;
T <- length(maturity)	
lambda<- 0.0609 # loading parameter

# Loading matrix 
source("Nelson.Siegel.factor.loadings.R")
Z <- Nelson.Siegel.factor.loadings(lambda=lambda,maturity=maturity)

# Betas and Yhat
source("Yhat.betas.R")
results <- Yhat.betas(Y=data) 
head(results$beta)
head(results$Yhat)

# VAR(1) coeffient matrix 
source("VARcoeff.R")
var<-VARcoeff(betas=results$beta) # Start point to pars$phi in Kalman-Filter-Dynamic-Nelson-Siegel

# It fits data 348. Obs.: I calculate the VAR coefficients matrix to the whole sample so that 
# I can use for forecasts in observation 349, for example.

betahat1 <- var[,2:4] %*% results$beta[347,]
Yhat1 <- Z %*% betahat1
ts.plot(t(data[348,])) # Observed
lines(Yhat1) # Fitted

# Compare VAR(1) with Package("vars") 
require(vars)
var1 <- VAR(results$beta, 1, type=c("const"),season = NULL, exogen = NULL)
var11<-summary(var1)
var11$varresult

# Compare with Package("YieldCurve")
require(YieldCurve)
maturity 	<- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
NSParameters<- Nelson.Siegel(rate=data, maturity=maturity)
head(NSParameters)
dim(NSParameters)
var2 <- VAR(NSParameters, 1, type=c("const"),season = NULL, exogen = NULL)
var22<-summary(var2)

# Comment:
# There is a difference between the betas of the Package("YieldCurve") and the betas of the function Yhat.betas
# since the Package has a "time-varying" loading parameter (lambda) for each yield curve observed, 
# and the function has only one fixed lambda (0.0609) for all observations.

# Proxy to time-varying volatility (h_{t}). See Koopman,Malle,Var der Wel(2010), p.342, Fig.4, Panel(B)
sigma<-c(rep(NA,348))
  for(i in 1:348)
    {
      sigma[i]<- var(as.numeric((data-results$Yhat)[i,]))
    }
ts.plot(sigma)

# Initials values for Kalman-Filter-Dynamic-Nelson-Siegel (Lyapunov equation) See Koopman,Malle,Var der Wel(2010), p.331
# Variance matrix
I <-diag(9)
K <-var[,2:4]
KK<-kronecker(K,K)
J <-var11$covres
JJ<-matrix(J,9,1)
V <-solve(I - KK) %*% JJ
V <-matrix(V,3,3)

# Start point to DNS-baseline with Kalman filter
para<-c(rep(NA,36))
# Start point to lambda 
para[1]<-lambda 
# Start point to pars$H 
para[2:18]<-sqrt(diag(var(data-results$Yhat))) 
# Start point to pars$phi
para[19]<-var[1,2]
para[20]<-var[1,3]
para[21]<-var[1,4]
para[22]<-var[2,2]
para[23]<-var[2,3]
para[24]<-var[2,4]
para[25]<-var[3,2]
para[26]<-var[3,3]
para[27]<-var[3,4]
# Start point to pars$mu 
para[28]<-mean(results$beta[,1]) 
para[29]<-mean(results$beta[,2]) 
para[30]<-mean(results$beta[,3])
# Start to pars$Q
QQ<-t(chol(var11$covres)) 
para[31]<-QQ[1,1]
para[32]<-QQ[2,1]
para[33]<-QQ[2,2]
para[34]<-QQ[3,1]
para[35]<-QQ[3,2]
para[36]<-QQ[3,3]
