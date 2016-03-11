
sink("Bayesian/NmixturePoissonRagged2m.jags")

cat("
    model {
    #Compute intensity for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Times){
    
    #Process Model
    log(lambda[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j] + beta2[i] * resources[i,j,k] + beta3[i] * Traitmatch[i,j] * resources[i,j,k]
    
    #True number of interactions
    N[i,j,k] ~ dpois(lambda[i,j,k])
    }
    }
    }

    #Observed counts for each day of sampling
    for (x in 1:Nobs){
    
    #Observation Process for cameras
    detect_cam[x]<-dcam[Bird[x]] * cam_surveys[x]

    #Observation Process for transects
    detect_transect[x]<-dtrans[Bird[x]] * trans_surveys[x]

    Yobs_camera[x] ~ dbin(detect_cam[x],N[Bird[x],Plant[x],Time[x]])    
    Yobs_transect[x] ~ dbin(detect_transect[x],N[Bird[x],Plant[x],Time[x]])    

    #Assess Model Fit

    #Fit discrepancy statistics
    #eval[x]<-detect[Bird[x]]*N[Bird[x],Plant[x],Camera[x]]
    #E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
    
    #ynew[x]~dbin(detect[Bird[x]],N[Bird[x],Plant[x],Camera[x]])
    #E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
    
    }
    

    #Species level priors
    
    for (i in 1:Birds){
    #Detect priors, logit transformed

    #For Cameras
    logit(dcam[i]) <- dcam_logit[i]
    dcam_logit[i] ~ dnorm(dprior_cam,tau_dcam)

    #For Transects
    logit(dtrans[i]) <- dtrans_logit[i]
    dtrans_logit[i] ~ dnorm(dprior_trans,tau_dtrans)

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
    
    #Detection group prior
    dprior_cam ~ dnorm(0,0.5)
    dprior_trans ~ dnorm(0,0.5)

    # Group intercept variance
    tau_alpha ~ dgamma(0.0001,0.0001)
    sigma_int<-pow(1/tau_alpha,2) 
    
    #Group effect detect camera
    tau_dcam ~ dunif(0,10)
    sigma_dcam<-pow(1/tau_dcam,.5)
    
    #Group effect detect camera
    tau_dtrans ~ dunif(0,10)
    sigma_dtrans<-pow(1/tau_dtrans,.5)

    #Group Effect of traits
    tau_beta1 ~ dgamma(0.0001,0.0001)
    sigma_slope1<-pow(1/tau_beta1,.5)
    
    #Group Effect of Resources
    tau_beta2 ~ dgamma(0.0001,0.0001)
    sigma_slope2<-pow(1/tau_beta2,.5)
    
    #Group Effect of Resources * Traits
    tau_beta3 ~ dgamma(0.0001,0.0001)
    sigma_slope3<-pow(1/tau_beta3,.5)
}
    ",fill=TRUE)

sink()
