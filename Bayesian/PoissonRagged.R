
sink("Bayesian/PoissonRagged.jags")

cat("
    model {
    for (x in 1:Nobs){
    
      # True state model for the only partially observed true state    
 log(lambda[Bird[x],Plant[x],Time[x]]) <- alpha[Bird[x]] + beta1[Bird[x]] * traitmatch[x] + beta2[Bird[x]] * resources[x] + beta3[Bird[x]] * resources[x] * traitmatch[x]      
      Yobs[x] ~ dpois(lambda[Bird[x],Plant[x],Time[x]] )    
    }
    
    for (i in 1:Birds){
      alpha[i] ~ dnorm(intercept,tau_alpha)
      beta1[i] ~ dnorm(gamma1,tau_beta1)    
      beta2[i] ~ dnorm(gamma2,tau_beta2)    
      beta3[i] ~ dnorm(gamma3,tau_beta3)    
    }
    
    #Hyperpriors
    #Slope grouping
    gamma1~dnorm(0.001,0.001)
    gamma2~dnorm(0.001,0.001)
    gamma3~dnorm(0.001,0.001)
    
    #Intercept grouping
    intercept~dnorm(0.001,0.001)
    
    # Group variance
    tau_alpha ~ dgamma(0.0001,0.0001)
    sigma_int<-pow(1/tau_alpha,0.5) #Derived Quantity
    
    #Slope
    tau_beta1 ~ dgamma(0.0001,0.0001)
    tau_beta2 ~ dgamma(0.0001,0.0001)
    tau_beta3 ~ dgamma(0.0001,0.0001)
    
    sigma_slope1<-pow(1/tau_beta1,0.5)
    sigma_slope2<-pow(1/tau_beta2,0.5)
    sigma_slope3<-pow(1/tau_beta3,0.5)
    
    }
    ",fill=TRUE)

sink()
