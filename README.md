# Evidence_Statistics
R code for evidence statistics 
```r

library(ggplot2)
library(tidyr)
library(dplyr)
library(VGAM)


#1.Function info_criterion

###1)
CAL_DIC<- function(N, NSIMS,tmu,null_nll,null_init,alt_nll, alt_init){
    #This function is to get delta info_criterion values between null hypothesis, using null_mle function,
    #and althernative hypothesis, using alt_mle function,
    #starting with initial values null_init and alt_init
    #N: sample size
    #tmu: true mean of the population, using normal distribution
    #null_nll: negative loglike function of null distribution
    #null_init: initial parameter to initiate optim
    #alt_nll: negative loglike function of alterative disttribution
    #alt_init: initial parameter to initiate optim
    
    mles.mat <-matrix(0, nrow=NSIMS, ncol=3); #Step 1: make matrix that stores dIC values
    for(i in 1:NSIMS){
        
        dataset<-rnorm(N,tmu,sd=15)
        
        # Step 2: Get Negative sum of loglike functions for null and alternative functions
        ith.null.results <- optim(par=null_init, fn=null_nll, method="BFGS", data=dataset)
        ith.alt.results <- optim(par=alt_init, fn=alt_nll, method="BFGS", data=dataset)
        
        
        # Step 3: Get information criterions
        k=length(null_init)-length(alt_init)
        dAIC <- 2*(ith.null.results$value-ith.alt.results$value)+2*k
        dAICc <- 2*(ith.null.results$value-ith.alt.results$value)+2*k*(k+1)/(N-k-1)
        dBIC <-2*(ith.null.results$value-ith.alt.results$value)+log(N)*k
        mles.mat[i,] <-c(dAIC,dAICc,dBIC)
    }
    colnames(mles.mat) <-c("dAIC","dAICc","dBIC")
    return(mles.mat)
}
#2) Fuctinon (A)4-1.
null_A4 <- function(parms,data){
    # parms:  parameters for normal MLE, mean=fixed, 100 and sd=parms[1]
    sig<-exp(parms[1])
    dat.len <- length(data)
    llvec <- rep(0,dat.len)
    for(i in 1:dat.len){
        llvec[i] <- dnorm(data[i],mean=100, sd=sig, log=TRUE)
    }
    tot.negloglike <- -sum(llvec)
    return(tot.negloglike)
}

#3) Function (A)4-2.
alt_A4 <- function(parms,data){
    # parms:  parameters for normal MLE, mean=parms[1] and sd=parms[2]
    mu=parms[1]
    sig=exp(parms[2])
    dat.len <- length(data)
    llvec <- rep(0,dat.len)
    for(i in 1:dat.len){
        llvec[i] <- dnorm(data[i],mean=mu, sd=sig, log=TRUE)
    }	
    tot.negloglike <- -sum(llvec)
    return(tot.negloglike)
}


#4) Function calculating delta IC sampling from Laplace distribution


CAL_DIC_Lplc<- function(N, NSIMS,tmu,tsig,null_nll,null_init,alt_nll, alt_init){
    #This function is to get delta info_criterion values between null hypothesis, using null_mle function,
    #and althernative hypothesis, using alt_mle function,
    #starting with initial values null_init and alt_init
    #N: sample size
    #tmu: true mean of the population, using laplace distribution
    #tsig: true sigma of the population, using laplace distribution
    #null_nll: negative loglike function of null distribution
    #null_init: initial parameter to initiate optim
    #alt_nll: negative loglike function of alterative disttribution
    #alt_init: initial parameter to initiate optim
    
    mles.mat <-matrix(0, nrow=NSIMS, ncol=3); #Step 1: make matrix that stores dIC values
    for(i in 1:NSIMS){
        
        dataset<-rlaplace(n=10,location=tmu,scale=sqrt(tsig^2/2)) #rlaplace using package VGAM
        
        # Step 2: Get Negative sum of loglike functions for null and alternative functions
        ith.null.results <- optim(par=null_init, fn=null_nll, method="BFGS", data=dataset)
        ith.alt.results <- optim(par=alt_init, fn=alt_nll, method="BFGS", data=dataset)
        
        
        # Step 3: Get information criterions
        k=length(null_init)-length(alt_init)
        dAIC <- 2*(ith.null.results$value-ith.alt.results$value)+2*k
        dAICc <- 2*(ith.null.results$value-ith.alt.results$value)+2*k*(k+1)/(N-k-1)
        dBIC <-2*(ith.null.results$value-ith.alt.results$value)+log(N)*k
        mles.mat[i,] <-c(dAIC,dAICc,dBIC)
    }
    colnames(mles.mat) <-c("dAIC","dAICc","dBIC")
    return(mles.mat)
}



##===================3. Solution to problem (A)-6
#The alternative model with mu=115 is true, therefore dataset is generated using true mean=115

ndata <-c(4:15); #Simluation is from N=4 to N=30
NSIMS <- 500; #Number of simulation is 10000
error.rate <- matrix(0, nrow= length(ndata), ncol=4); #misleading information frequency, for AIC, AICc, BIC
c1=c(5); #intial values for nll mle, #parameter=1
c2=c(200,5); #intial values for alt mle, #parameter=2

for (N in min(ndata):max(ndata)){
    
    DIC=CAL_DIC(N,NSIMS,tmu=115,null_A4,c1,alt_A4,c2); #true mean is 115 in this dataset, alternative is the true model
    #Similate under the null, chance of your decision is wrong under the null model
    error.rate[N-3,] <- c(N,sum(DIC[,1]<2),sum(DIC[,2]<2),sum(DIC[,3]<2)); #K=2, frequency that exceeds 2 is recorded or given N=# of samples
    # Step 4: plot a histogram of the mles and compare to true values
}
colnames(error.rate)=c("N","AIC","AICc","SIC") #columns are labeled
error.rate=as.data.frame(error.rate)
ans.A.6<-error.rate %>% gather("N") #Changing the matrix into dataframe for the final answer
ans.A.6$value=ans.A.6$value/NSIMS #Convert frequency to probability
colnames(ans.A.6)=c("N","dIC","prob") #label each columes
#Plot the final answer
ggplot(ans.A.6, aes(N,prob))+geom_line(aes(colour=dIC))+ylab("Misleading Evidence")+xlab("Number of Samples")+ggtitle("Simulations of different ICs (10,000 times)")





```
