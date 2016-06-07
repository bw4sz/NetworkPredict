
sink("Bayesian/NmixturePoissonRagged2m.jags")

cat("
    model {
    #Compute true state for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Times){
    
    #Process Model
    logit(rho[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j] + beta2[i] * resources[i,j,k]
    
    #True State
    S[i,j,k] ~ dbern(rho[i,j,k])
    }
    }
    }
    
    #Observation Model
    for (x in 1:Nobs){
    
    #Observation Process for cameras
    detect_cam[x]<-dcam * cam_surveys[x]

    #Observation Process for transects
    detect_transect[x]<-dtrans * trans_surveys[x]

    Yobs_camera[x] ~ dbern(detect_cam[x] * S[Bird[x],Plant[x],Time[x]])    
    Yobs_transect[x] ~ dbern(detect_transect[x] * S[Bird[x],Plant[x],Time[x]])    

    #Assess Model Fit

    #Fit discrepancy statistics
    #eval[x]<-detect[Bird[x]]*S[Bird[x],Plant[x],Camera[x]]
    #E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
    
    #ynew[x]~dbin(detect[Bird[x]],S[Bird[x],Plant[x],Camera[x]])
    #E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    
    }

    #Transect Prior
    #Detect priors, logit transformed
    #For Cameras
    logit(dcam) <- dcam_logit
    dcam_logit ~ dnorm(0,0.386)
    
    #For Transects
    logit(dtrans) <- dtrans_logit
    dtrans_logit ~ dnorm(0,0.386)
    
    #Species level priors
    for (i in 1:Birds){
      alpha[i] ~ dnorm(intercept,tau_alpha)
      beta1[i] ~ dnorm(gamma1,tau_beta1)    
      beta2[i] ~ dnorm(gamma2,tau_beta2)    
    }

    #Hyperpriors
    #Slope grouping
    gamma1~dnorm(0,0.0001)
    gamma2~dnorm(0,0.0001)
    
    #Intercept grouping
    intercept~dnorm(0,0.0001)
    
    #Detection group prior
    #dprior_cam ~ dnorm(0,0.386)
    #dprior_trans ~ dnorm(0,0.386)

    # Group intercept variance
    tau_alpha ~ dt(0,1,1)I(0,1)
    sigma_alpha<-pow(1/tau_alpha,2) 
    
    #Group effect detect camera
    #tau_dcam ~ dunif(0,10)
    #sigma_dcam<-pow(1/tau_dcam,.5)
    
    #Group effect detect camera
    #tau_dtrans ~ dunif(0,10)
    #sigma_dtrans<-pow(1/tau_dtrans,.5)

    #Group Effect of traits
    tau_beta1 ~ dt(0,1,1)I(0,)
    sigma_slope1<-pow(1/tau_beta1,.5)
    
    #Group Effect of Resources
    tau_beta2 ~ dt(0,1,1)I(0,)
    sigma_slope2<-pow(1/tau_beta2,.5)

}
    ",fill=TRUE)

sink()
