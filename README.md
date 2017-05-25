# Evidence_Statistics
R code for evidence statistics 
```r

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





```
