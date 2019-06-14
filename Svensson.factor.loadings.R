# Loading matrix Lambda function
Svensson.factor.loadings<-function(lambda1,lambda2,maturity) 
	{
	  column1<-rep.int(1,length(maturity))
	  column2<-(1-exp(-lambda1*maturity))/(lambda1*maturity)
	  column3<-column2-exp(-lambda1*maturity)
	  column4<-(1-exp(-lambda2*maturity))/(lambda2*maturity)-exp(-lambda2*maturity)
	  
	  lambmat<-cbind(column1,column2,column3,column4)
	  
	  lambmat
	}  
