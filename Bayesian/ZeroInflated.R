
sink("ZeroInflated.jags")

cat("
    model {
    for (i in 1:Birds){
      for (j in 1:Plants){
        for (k in 1:Months){

      # True state model for the only partially observed true state    
      log(lambda[i,j,k])<- alpha[i] + beta1[i] * traitmatch[i,j,k] + beta2[i] * resources[i,j,k] + beta3[i] * resources[i,j,k] * traitmatch[i,j,k]
      N[i,j,k] ~ dpois(lambda[i,j,k])

      # Observation model for the actual observations
      p~ dbin(detect[i]) 
      Y[i,j,k] <- N[i,j,k] * p 
      
      #Fit discrepancy statistics
      eval[i,j,k]<-detect[i]*N[i,j,k]
      E[i,j,k]<-pow((Y[i,j,k]-eval[i,j,k]),2)/(eval[i,j,k]+0.5)

      y.new[i,j,k]~dbin(detect[i],N[i,j,k])
      E.new[i,j,k]<-pow((y.new[i,j,k]-eval[i,j,k]),2)/(eval[i,j,k]+0.5)
        }
      }
    }
    
    for (i in 1:Birds){
    detect[i] ~ dunif(0,1) # Detection for each bird species
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
    
    
    #derived posterior check
    fit<-sum(E[,,]) #Discrepancy for the observed data
    fitnew<-sum(E.new[,,])
    
    }
    ",fill=TRUE)

sink()
