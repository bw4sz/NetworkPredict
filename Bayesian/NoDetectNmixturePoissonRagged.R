
sink("Bayesian/NoDetectNmixturePoissonRagged.jags")

cat("
    model {

    for (i in 1:Birds){
      for (j in 1:Plants){
        for (k in 1:Times){
          #Process Model
              log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j] 
                gamma[i,j,k]=beta2[i] * resources[i,j,k] 
        }
      }
    }
    
    for (x in 1:Nobs){
    
    # Covariates for observed state   
    Yobs[x] ~ dpois(lambda[Bird[x],Plant[x],Time[x]] * gamma[i,j,k])    
    
    #Assess Model Fit
    
    #Fit discrepancy statistics
    eval[x]<-lambda[Bird[x],Plant[x],Time[x]]
    E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
    
    ynew[x]~dpois(lambda[Bird[x],Plant[x],Time[x]])
    E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    
    }
    
    for (i in 1:Birds){
    alpha[i] ~ dnorm(alpha_mu,alpha_tau)
    beta1[i] ~ dnorm(beta1_mu,beta1_tau)
    beta2[i] ~ dnorm(beta2_mu,beta2_tau)    
    }
    
    #Hyperpriors
    #Slope grouping
    beta1_mu~dnorm(0,0.0001)
    beta2_mu~dnorm(0,0.0001)

    #Intercept grouping
    alpha_mu~dnorm(0,0.0001)

    # Group intercept variance
    alpha_tau ~ dgamma(0.0001,0.0001)
    alpha_sigma<-pow(1/alpha_tau,0.5) 
    
    #Derived Quantity
    
    #Slope variance, turning precision to sd
    beta1_tau ~ dgamma(0.0001,0.0001)
    beta1_sigma<-pow(1/beta1_tau,0.5)
    
    beta1_tau ~ dgamma(0.0001,0.0001)
    beta1_sigma<-pow(1/beta1_tau,0.5)

    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()
