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
#setwd("C:\\...")
require(xts)
#data <-read.csv("bonds.csv",header = TRUE, sep = ";")
data <- as.matrix(data[,which(names(data)=="M3"):which(names(data)=="M120")])
datas <- seq(as.Date("1972/1/1"), by = "month", length.out = 348)
data  <- xts(data, order.by = datas)
head(data)
dim(data)

l <-dim(data)[1] # observations
T <- length(maturity)	

maturity<-c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120) # 17 maturities in the data;
lambda<- 0.0609 # loading parameter

# Loading matrix Lambda function
Nelson.Siegel.factor.loadings<-function(lambda,maturity) 
	{
	  column1<-rep.int(1,length(maturity))
	  column2<-(1-exp(-lambda*maturity))/(lambda*maturity)
	  column3<-column2-exp(-lambda*maturity)
	  
	  lambmat<-cbind(column1,column2,column3)
	  
	  lambmat
	}  

Z <- Nelson.Siegel.factor.loadings(lambda,maturity) # loading matrix

# Betas and Yhat function
NS <- function(Y)
	{
	beta <- solve(t(Z) %*% Z) %*% t(Z) %*% t(Y) # Z=17x3, Z'=3x17, Y=348x17, Y'=17x348
	beta <- t(beta)
	Yhat <- beta %*% t(Z)

	return(list(Yhat=Yhat,beta=beta))
	}

results <- NS(Y=data) 

head(results$beta) # betas
head(results$Yhat) # Yhat

# Start to pars$mu in Kalman-Filter-Dynamic-Nelson-Siegel
mean(results$beta[,1]) 
mean(results$beta[,2]) 
mean(results$beta[,3])

# Residuals
res <- data-results$Yhat
head(res)

# Start to pars$H in Kalman-Filter-Dynamic-Nelson-Siegel (variance matrix of residuals)
sqrt(diag(var(res))) 

# Proxy to time-varying volatility (h_{t}). See Koopman,Malle,Var der Wel(2010), p.342, Fig.4, Panel(B)
sigma<-c(rep(NA,348))
for(i in 1:348)
{
sigma[i]<- var(as.numeric(res[i,]))
}
ts.plot(sigma)

# Vector autoregressive coeffient matrix VAR(1)
YY <-results$beta
XX <- cbind(rep(1,length(results$beta[,1])),YY)
XX <-XX[-l,] # exclude observation 348
YY <-YY[-1,] # exclude observation 1
var <- solve(t(XX) %*% XX) %*% t(XX) %*% YY
var <-t(var)
var # Start to pars$phi in Kalman-Filter-Dynamic-Nelson-Siegel

# Fitting the observation 348. Obs.: The VAR calculation was done on the whole sample,
# so it would be a forecast on yields in observation 349.
betahat1 <- var[,2:4] %*% results$beta[347,]
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
# Variance matrix of betas
I <-diag(9)
K <-var[,2:4]
KK<-kronecker(K,K)
J <-var11$covres
JJ<-matrix(J,9,1)
V <-solve(I - KK) %*% JJ
V <-matrix(V,3,3)

# Compare with Package("YieldCurve")
require(YieldCurve)
maturity 	<- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
NSParameters<- Nelson.Siegel( rate=data, maturity=maturity)
head(NSParameters)
dim(NSParameters)
var2 <- VAR(NSParameters, 1, type=c("const"),season = NULL, exogen = NULL)
var22<-summary(var2)

# Comment:
# There is a difference between the betas of the Package("YieldCurve") and the betas of the function 
# since the Package has a "time-varying" loading parameter (lambda) for each observed yield curve, 
# and the function has only one fixed lambda (0.0609) for the observations.
