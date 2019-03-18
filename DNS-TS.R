# Autor: Werley Cordeiro

list.of.packages <- c("highfrequency","vars","YieldCurve")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(YieldCurve)
library(vars)
library(highfrequency)



data <-read.csv("https://www.dropbox.com/s/inpnlugzkddp42q/bonds.csv?dl=1",header = TRUE, sep = ";")
#setwd("C:\\Users\\werle\\Dropbox\\Mestrado_UFSC\\4_Semestre\\DNS-TS")
require(xts)
#data <-read.csv("bonds.csv",header = TRUE, sep = ";")
data <- as.matrix(data[,which(names(data)=="M3"):which(names(data)=="M120")])
datas <- seq(as.Date("1972/1/1"), by = "month", length.out = 348)
data  <- xts(data, order.by = datas)


l <-dim(data)[1] 

maturity<-c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120) # Primeiro: definir as maturidades que tenho nos dados;
lambda<- 0.0609 # Segundo: definir o valor do parâmetro de decaimento;
T <- length(maturity)	

	Nelson.Siegel.factor.loadings<-function(lambda,maturity)
	{
	  column1<-rep.int(1,length(maturity))
	  column2<-(1-exp(-lambda*maturity))/(lambda*maturity)
	  column3<-column2-exp(-lambda*maturity)
	  
	  lambmat<-cbind(column1,column2,column3)
	  
	  lambmat
	}  

Z <- Nelson.Siegel.factor.loadings(lambda,maturity)	# Terceiro: calcular matriz dos fatores de carregamento;

NS <- function(Y){

beta <- solve(t(Z) %*% Z) %*% t(Z) %*% t(Y) # Z=17x3, Z'=3x17, Y=348x17, Y'=17x348
beta <- t(beta)
Yhat <- beta %*% t(Z)
	
	return(list(Yhat=Yhat,beta=beta))
}

# Quarto: calcular os betas OLS das taxas;
results <- NS(Y=data) 
head(results$beta)
dim(results$beta)

head(results$Yhat)
dim(results$Yhat)

# Start para pars$mu no DNS-baseline (FK)
mean(results$beta[,1]) 
mean(results$beta[,2]) 
mean(results$beta[,3])

# Resíduo
res <- data-results$Yhat
head(res)
dim(res)

sqrt(diag(var(res))) # Start para pars$H no DNS-baseline (FK) (matriz variâncias diagonal)

# Proxy to time-varying volatility (h_{t})
sigma<-c(rep(NA,348))
for(i in 1:348)
{
sigma[i]<- var(as.numeric(res[i,]))
}
ts.plot(sigma)

# Estimates of latent factors VAR(1) model
YY <-results$beta
XX <- cbind(rep(1,length(results$beta[,1])),YY)
XX <-XX[-l,]
YY <-YY[-1,]
var <- solve(t(XX) %*% XX) %*% t(XX) %*% YY
var <-t(var)
var #Start para pars$phi 

# Forecast one step ahead Obs.: O cálculo do VAR foi feito sobre toda a amostra, portanto é uma "previsão"...
betahat1 <- var[,2:4] %*% results$beta[347,]
Yhat1 <- Z %*% betahat1
ts.plot(t(data[348,])) # Observed
lines(Yhat1) # Forecast

# Compare Package("vars") 
require(vars)
var1 <- VAR(results$beta, 1, type=c("const"),season = NULL, exogen = NULL)
var11<-summary(var1)
t(chol(var11$covres)) #Start para pars$Q ******************************

# Initials values for DNS-baseline (Lyapunov equation)
I <-diag(9)
K <-var[,2:4]
KK<-kronecker(K,K)
J <-var11$covres
JJ<-matrix(J,9,1)
V <-solve(I - KK) %*% JJ
V <-matrix(V,3,3)

# Compare Package("YieldCurve")
require(YieldCurve)
maturity 	<- c(3,6,9,12,15,18,21,24,30,36,48,60,72,84,96,108,120)
NSParameters<- Nelson.Siegel( rate=data, maturity=maturity)
NSParameters
var2 <- VAR(NSParameters, 1, type=c("const"),season = NULL, exogen = NULL)
var22<-summary(var2)
head(var22)
dim(var22)

# Comentário:
# Há diferença entre os betas do pacote (YieldCurve) e os betas da função, pois o pacote tem um lambda para cada ETTJ, e a função 
# tem apenas um lambda fixo (0.0609) para todas as ETTJ's.



