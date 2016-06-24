
sink("Bayesian/NmixturePoissonRagged2m.jags")

cat("
    model {
    #Compute true state for each pair of birds and plants
    for (i in 1:Birds){
    for (j in 1:Plants){
    for (k in 1:Times){
    
    #Process Model
    logit(rho[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j] + beta2[i] * resources[i,j,k] + beta3[i] * resources[i,j,k] * Traitmatch[i,j]
    
    #True State
    S[i,j,k] ~ dbern(rho[i,j,k])
    }
    }
    }
    
    #Observation Model
    for (x in 1:Nobs){
    
    #Observation Process for cameras
    detect_cam[x]<-dcam[Bird[x]] * cam_surveys[x]
    
    #Observation Process for transects
    detect_transect[x]<-dtrans[Bird[x]] * trans_surveys[x]
    
    Yobs_camera[x] ~ dbern(detect_cam[x] * S[Bird[x],Plant[x],Time[x]])    
    Yobs_transect[x] ~ dbern(detect_transect[x] * S[Bird[x],Plant[x],Time[x]])    
    
    #     #Assess Model Fit - Posterior predictive check
    # 
    #     #Fit discrepancy statistics
    #     #Camera
    #     eval_cam[x]<-detect_cam[Bird[x]]*S[Bird[x],Plant[x],Time[x]]
    #     E_cam[x]<-pow((Yobs_camera[x]-eval_cam[x]),2)/(eval[x]+0.5)
    #     
    #     ynew_cam[x]~dbern(detect_cam[Bird[x]] * S[Bird[x],Plant[x],Time[x]])
    #     E.new_cam[x]<-pow((ynew_cam[x]-eval_cam[x]),2)/(eval_cam[x]+0.5)
    # 
    #     #Transect
    #     eval_transect[x]<-detect_transect[Bird[x]]*S[Bird[x],Plant[x],Time[x]]
    #     E_transect[x]<-pow((Yobs_transect[x]-eval_transect[x]),2)/(eval[x]+0.5)
    #     
    #     ynew_transect[x]~dbern(detect_transect[Bird[x]] * S[Bird[x],Plant[x],Time[x]])
    #     E.new_transect[x]<-pow((ynew_transect[x]-eval_transect[x]),2)/(eval_transect[x]+0.5)
    #     
    #     fit_trans<-sum(E_transect[]) #Discrepancy for the observed data
    #     fitnew_trans<-sum(E.new_transect[])
    # 
    #     fit_cam<-sum(E_cam[]) #Discrepancy for the observed data
    #     fitnew_cam<-sum(E.new_cam[])
    
    }
    
    #Priors
    #Observation model
    #Detect priors, logit transformed - Following lunn 2012 p85
    
    for(x in 1:Birds){
    #For Cameras
    dcam[x] ~ dnorm(dprior_cam,tau_dcam)
    
    #For Transects
    dtrans[x] ~ dnorm(dprior_trans,tau_dtrans)
    }
    
    #Detection group prior
    dprior_cam ~ dnorm(0,0.386)
    dprior_trans ~ dnorm(0,0.386)
    
    #Group effect detect camera
    tau_dcam ~ dunif(0,1000)
    sigma_dcam<-pow(1/tau_dcam,.5)
    
    #Group effect detect camera
    tau_dtrans ~ dunif(0,1000)
    sigma_dtrans<-pow(1/tau_dtrans,.5)
    
    #Process Model
    #Species level priors
    for (i in 1:Birds){
    
    #Intercept
    #logit prior, then transform for plotting
    alpha[i] ~ dnorm(alpha_mu,alpha_tau)
    
    #Traits slope 
    beta1[i] ~ dnorm(beta1_mu,beta1_tau)    
    
    #Plant slope
    beta2[i] ~ dnorm(beta2_mu,beta2_tau)    
    
    #Interaction slope
    beta3[i] ~ dnorm(beta3_mu,beta3_tau)    
    }
    
    #Group process priors
    
    #Intercept 
    alpha_mu ~ dnorm(0,0.386)
    alpha_tau ~ dunif(0,1000)
    alpha_sigma<-pow(1/alpha_tau,0.5) 
    
    #Trait
    beta1_mu~dnorm(0,0.386)
    beta1_tau ~ dunif(0,1000)
    beta1_sigma<-pow(1/beta1_tau,0.5)
    
    #Resources
    beta2_mu~dnorm(0,0.386)
    beta2_tau ~ dunif(0,1000)
    beta2_sigma<-pow(1/beta2_tau,0.5)
    
    #Interaction
    beta3_mu~dnorm(0,0.386)
    beta3_tau ~ dunif(0,1000)
    beta3_sigma<-pow(1/beta3_tau,0.5)
    
    }
    ",fill=TRUE)

sink()