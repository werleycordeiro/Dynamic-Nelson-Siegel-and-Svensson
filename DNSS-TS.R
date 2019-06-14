# author: Werley Cordeiro
# werleycordeiro@gmail.com

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
lambda1<- 0.0609 # loading parameter 1
lambda2<- 0.1218 # loading parameter 2

# Loading matrix 
source("Svensson.factor.loadings.R")
Z <- Svensson.factor.loadings(lambda1=lambda1,lambda2=lambda2,maturity=maturity)

# Betas and Yhat function
source("Yhat.betas.R")
results <- NS(Y=data) 
head(results$beta)
head(results$Yhat)

# Start point to pars$mu 
mean(results$beta[,1]) 
mean(results$beta[,2]) 
mean(results$beta[,3])
mean(results$beta[,4])

# Start point to pars$H 
sqrt(diag(var(data-results$Yhat))) 

# Proxy to time-varying volatility (h_{t}). See Koopman,Malle,Var der Wel(2010), p.342, Fig.4, Panel(B)
sigma<-c(rep(NA,348))
  for(i in 1:348)
    {
      sigma[i]<- var(as.numeric((data-results$Yhat)[i,]))
    }
ts.plot(sigma)

# VAR(1) coeffient matrix 
source("VARcoeff.R")
var<-VARcoeff(betas=results$beta) # Start point to pars$phi 

# Fitting the observation 348. Obs.: The VAR calculation was done on the whole sample,
# so it would be used for forecasts on yields in observation 349.

betahat1 <- var[,2:5] %*% results$beta[347,]
Yhat1 <- Z %*% betahat1
ts.plot(t(data[348,])) # Observed
lines(Yhat1) # Fitted

# Compare VAR(1) with Package("vars") 
require(vars)
var1 <- VAR(results$beta, 1, type=c("const"),season = NULL, exogen = NULL)
var11<-summary(var1)
var11$varresult

# Start to pars$Q in Kalman-Filter-Dynamic-Nelson-Siegel
t(chol(var11$covres)) 

# Initials values for Kalman-Filter-Dynamic-Nelson-Siegel (Lyapunov equation) See Koopman,Malle,Var der Wel(2010), p.331
# Variance matrix
I <-diag(16)
K <-var[,2:5]
KK<-kronecker(K,K)
J <-var11$covres
JJ<-matrix(J,16,1)
V <-solve(I - KK) %*% JJ
V <-matrix(V,4,4)

# Compare with Package("YieldCurve")
require(YieldCurve)
maturity 	<- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
NSParameters<- Svensson(rate=data, maturity=maturity)
head(NSParameters)
dim(NSParameters)
var2 <- VAR(NSParameters, 1, type=c("const"),season = NULL, exogen = NULL)
var22<-summary(var2)

# Comment:
# There is a difference between the betas of the Package("YieldCurve") and the betas of the function Yhat.betas
# since the Package has a "time-varying" loading parameter (lambda) for each yield curve observed, 
# and the function has only one fixed lambda (0.0609) for all observations.