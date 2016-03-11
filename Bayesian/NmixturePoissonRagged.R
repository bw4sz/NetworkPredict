
sink("Bayesian/NmixturePoissonRagged.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Cameras){
    
    #Process Model
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j] + beta2[i] * resources[i,k] + beta3[i] * Traitmatch[i,j] * resources[i,k]
    
    #For each camera - there is a latent count
    N[i,j,k] ~ dpois(lambda[i,j,k])
    }
    }
    }
    
    
    #Observed counts for each day of sampling at that camera
    for (x in 1:Nobs){
    
    #Observation Process
    Yobs[x] ~ dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])    
    
    #Assess Model Fit
    
    #Fit discrepancy statistics
    eval[x]<-detect[Bird[x]]*N[Bird[x],Plant[x],Camera[x]]
    E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
    
    ynew[x]~dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])
    E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    
    }
    

    #Species level priors

    for (i in 1:Birds){
    detect[i] ~ dunif(0,0.5)
    alpha[i] ~ dnorm(intercept,tau_alpha)
    beta1[i] ~ dnorm(gamma1,tau_beta1)    
    beta2[i] ~ dnorm(gamma2,tau_beta2)    
    beta3[i] ~ dnorm(gamma3,tau_beta3)    
    }
    
    #Hyperpriors
    #Slope grouping
    gamma1~dnorm(0,0.0001)
    gamma2~dnorm(0,0.0001)
    gamma3~dnorm(0,0.0001)

    
    #Intercept grouping
    intercept~dnorm(0,0.0001)
    
    # Group intercept variance
    tau_alpha ~ dgamma(0.0001,0.0001)
    sigma_int<-pow(1/tau_alpha,0.5) 
    
    #Derived Quantity
    
    #Slope variance, turning precision to sd

    #Group Effect of traits
    tau_beta1 ~ dgamma(0.0001,0.0001)
    sigma_slope1<-pow(1/tau_beta1,0.5)

    #Group Effect of Resources
    tau_beta2 ~ dgamma(0.0001,0.0001)
    sigma_slope2<-pow(1/tau_beta2,0.5)

    #Group Effect of Resources * Traits
    tau_beta3 ~ dgamma(0.0001,0.0001)
    sigma_slope3<-pow(1/tau_beta3,0.5)
    
    #derived posterior check
    fit<-sum(E[]) #Discrepancy for the observed data
    fitnew<-sum(E.new[])
    
    }
    ",fill=TRUE)

sink()
