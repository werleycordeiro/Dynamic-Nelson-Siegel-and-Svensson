# Loading matrix Lambda function
Nelson.Siegel.factor.loadings<-function(lambda,maturity) 
	{
	  column1<-rep.int(1,length(maturity))
	  column2<-(1-exp(-lambda*maturity))/(lambda*maturity)
	  column3<-column2-exp(-lambda*maturity)
	  
	  lambmat<-cbind(column1,column2,column3)
	  
	  lambmat
	}  
