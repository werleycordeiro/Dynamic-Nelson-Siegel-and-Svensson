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
results <- Yhat.betas(Y=data) 
head(results$beta)
head(results$Yhat)

# VAR(1) coeffient matrix 
source("VARcoeff.R")
var<-VARcoeff(betas=results$beta) # Start point to pars$phi 

# It fits data 348. Obs.: I calculate the VAR coefficients matrix to the whole sample so that 
# I can use for forecasts in observation 349, for example.

betahat1 <- var[,2:5] %*% results$beta[347,]
Yhat1 <- Z %*% betahat1
ts.plot(t(data[348,])) # Observed
lines(Yhat1) # Fitted

# Compare VAR(1) with Package("vars") 
require(vars)
var1 <- VAR(results$beta, 1, type=c("const"),season = NULL, exogen = NULL)
var11<-summary(var1)
var11$varresult

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

# Proxy to time-varying volatility (h_{t}). See Koopman,Malle,Var der Wel(2010), p.342, Fig.4, Panel(B)
sigma<-c(rep(NA,348))
  for(i in 1:348)
    {
      sigma[i]<- var(as.numeric((data-results$Yhat)[i,]))
    }
ts.plot(sigma)

# Comment:
# There is a difference between the betas of the Package("YieldCurve") and the betas of the function Yhat.betas
# since the Package has a "time-varying" loading parameter (lambda) for each yield curve observed, 
# and the function has only one fixed lambda (0.0609 and 0.1218) for all observations.

# Start point to DNSS-baseline with Kalman filter
para<-c(rep(NA,49))
# Start point to lambda 
para[1]<-lambda1
para[2]<-lambda2
# Start point to pars$H 
para[3:19]<-sqrt(diag(var(data-results$Yhat))) 
# Start point to pars$phi
para[20]<-var[1,2]
para[21]<-var[1,3]
para[22]<-var[1,4]
para[23]<-var[1,5]
para[24]<-var[2,2]
para[25]<-var[2,3]
para[26]<-var[2,4]
para[27]<-var[2,5]
para[28]<-var[3,2]
para[29]<-var[3,3]
para[30]<-var[3,4]
para[31]<-var[3,5]
para[32]<-var[4,2]
para[33]<-var[4,3]
para[34]<-var[4,4]
para[35]<-var[4,5]
# Start point to pars$mu 
para[36]<-mean(results$beta[,1]) 
para[37]<-mean(results$beta[,2]) 
para[38]<-mean(results$beta[,3])
para[39]<-mean(results$beta[,4])
# Start to pars$Q
QQ<-t(chol(var11$covres)) 
para[40]<-QQ[1,1]
para[41]<-QQ[2,1]
para[42]<-QQ[2,2]
para[43]<-QQ[3,1]
para[44]<-QQ[3,2]
para[45]<-QQ[3,3]
para[46]<-QQ[4,1]
para[47]<-QQ[4,2]
para[48]<-QQ[4,3]
para[49]<-QQ[4,4]
