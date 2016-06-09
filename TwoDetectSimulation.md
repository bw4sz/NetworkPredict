# Simulated Data for Two Detection Methods for Observing Species Interactions
Ben Weinstein  
December 25, 2015  

#Summary





```r
#load("Simulation_2M.RData")
```

There is some underlying network of interactions between hummingbird species  i and plant species j. We observe these interactions using transects across elevation ranges and cameras at individual flowers. To combine these data to jointly estimate the importance of trait-matching and resources on interaction intensity we need a hierarchical occupancy model that accounts for 1) the difference in sampling effort between survey types, 2) The variable number of replicates per species, 3) The difference in detectability of interactions based on survey type. The occupancy model below uses months as our estimated latent state. There are two surveys per month, and a variable number of cameras for each flower, often with no cameras on a given flower in a month.

# True Simulated values

Hummingbird Species =5

Plant Species=6

Survey Periods = 50

Detection Probability for Camera = 0.25

Detection Probability for Transect = 0.6

Hummingbird species are connected by the following hyperpriors

* intercept<-1

* alpha_sigma<- 0.2

Effect of Trait-matching

* gamma1=-0.5

* sigma_slope1<- 0.2

Effect of Resources

* gamma2=0

* sigma_slope2<- 0.1

Interaction effect of resources * traitmatch

* gamma3=0.3

* sigma_slope3<- 0.1

Bill sizes
Bill<-rpois(h_species,10)

Corolla sizes

Corolla<-rpois(plant_species,15)

Survey periods are 70% cameras, 30% Transect
Transects have two replicates.Cameras have variable number of replicates, modeled as rpois(lambda=0.5).

Resources are scored as either 'High' or 'Low' and is modeled as rbinom(n=1,size=1,prob=0.5)


```r
h_species=5
plant_species=6
Times=50
detection_cam=0.25
detection_trans=0.6

#which records are camera, which are transects?
mt<-rbinom(Times,1,0.7)
mt[which(mt==1)]<-"Camera"
mt[!mt=="Camera"]<-"Transect"

#Bill sizes
Bill<-rpois(h_species,10)

#Corolla sizes
Corolla<-rpois(plant_species,15)

#Subtract both and take absolute value
traitmatch<-abs(sapply(Corolla,function(x) x - Bill))

#fill out for each month
traitarray<-array(NA,dim=c(h_species,plant_species,Times))
#fill for each month
for (x in 1:Times){
  traitarray[,,x]<-traitmatch 
}

#simulate some poisson distributed resource counts for each replicate
#this will be same for each species to start with.
resources<-array(NA,dim=c(h_species,plant_species,Times))

#fill for each month
for (x in 1:Times){
  resources[,,x]<-rpois(1,10)  
}

#standardize predictors
#resources<-array(data=scale(resources,center=TRUE,scale=TRUE),dim=c(h_species,plant_species,Times))

#regression slope for trait-matching and resources

#Intercept
alpha_mu<-2
alpha_sigma<- 0.05

#trait match
beta1_mu=-0.5
beta1_sigma<- 0.05

#resources
beta2_mu=0
beta2_sigma<- 0.05

#loop through each species and plants

#draw values from hierarcichal distributions
beta1<-rnorm(h_species,beta1_mu,beta1_sigma)
beta2<-rnorm(h_species,beta2_mu, beta2_sigma)

alpha<-rnorm(h_species,alpha_mu,alpha_sigma)

phi<-inv.logit(alpha + beta1 * traitarray + beta2 * resources)

#How many cameras for each flower during each time period?
true_interactions<-array(data=sapply(phi,function(x){rbinom(1,1,prob=x)}),dim=c(h_species,plant_species,Times))

#combine and melt into a single datafFrame
mdat<-dcast(melt(list(y=true_interactions,traitmatch=traitarray,resources=resources)),Var1+Var2+Var3~L1)
colnames(mdat)<-c("Bird","Plant","Time","resources","traitmatch","True_state")

#Merge the survey type
mdat<-merge(mdat,data.frame(Time=1:Times,Survey_Type=mt))

##Observation models
dat<-list()
  
for (x in 1:nrow(mdat)){
  if(mdat$Survey_Type[x]=="Transect"){
    df<-data.frame(Y_Transect=rbinom(2,mdat$True_state[x],prob=detection_trans))
    dat[[x]]<-cbind(mdat[x,],df)
  } else{
        cams<-rpois(1,0.4)
        if(cams==0){next}
        df<-data.frame(Y_Camera=rbinom(cams,mdat$True_state[x],prob=detection_cam))
        dat[[x]]<-cbind(mdat[x,],df)
  }
}

mdat<-rbind_all(dat)
```

# Observed Data


```r
mdatm<-melt(mdat,measure.vars = c("True_state","Y_Camera","Y_Transect"))

ggplot(mdatm,aes(x=traitmatch,y=value,col=variable)) + geom_point(alpha=.5) + geom_smooth(method="glm",method.args=list(family="binomial"),linetype="dashed",size=1.1) + ggtitle("Correlation in Simulated Data") + labs(x="Difference in Bill and Corolla Length (mm)",y="Probability of Interactions",col="Observation Process") + theme_bw()
```

![](TwoDetectSimulation_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
#traitmatch dataframe
Traitmatch<-mdat %>% group_by(Bird,Plant) %>% summarize(v=unique(traitmatch)) %>% acast(Bird~Plant,value.var="v")

TimeResources<-mdat %>% group_by(Bird,Time,Plant) %>% summarize(v=unique(resources)) %>% acast(Bird~Plant~Time,value.var="v",fill=0)
```

#Hierarchical Occupancy Model

For hummingbird species i feeding on plant species j observed at time k and sampling event d. 

Observation Model

Transects

$$ Y_{Transect_{i,j,k,d}} \sim Bernoulli(S_{i,j,k} * \omega_{Transect})$$
$$ \omega_{Transect} <- \phi_{Transect}* EffortTransect_k $$

Cameras

$$ Y_{Camera_{i,j,k,d}} \sim Binomial(S_{i,j,k} * \omega_{Camera})$$
$$ \omega_{Camera} <- \phi_{Camera} * EffortCamera_k $$

Process Model

$$ S_{i,j,k} \sim Binomial(\rho_{i,j,k}) $$
$$ logit(\rho_{i,j,k}) = \alpha_i + \beta_{1,i} * Traitmatch_{i,j} + \beta_{2,i} *Resources_{j,k} $$


**Priors**

$$ \phi_{Camera} \sim Uniform(0,1) $$
$$ \phi_{Transect} \sim Uniform(0,1) $$
$$\alpha_i \sim Normal(\mu_\alpha,\tau_{\alpha})$$
$$\beta_{1,i} \sim Normal(\mu_{\beta_1},\tau_{\beta_1})$$
$$\beta_{2,i} \sim Normal(\mu_{\beta_2},\tau_{\beta_2})$$

**Hyperpriors**

Group Level Means

$$\mu_{\beta_1} \sim Normal(0,1.67)$$
$$\mu_{\beta_2} \sim Normal(0,1.67)$$
$$ \mu_{\alpha} \sim Normal(0,1.67)$$

Group Level Variance

$$\tau_{\alpha} \sim Uniform(0,100)$$
$$\tau_{\beta_1} \sim Uniform(0,100)$$
$$\tau_{\beta_2} \sim Uniform(0,100)$$

**Derived quantities**

$$\sigma_{\alpha} = \sqrt[2]{\frac{1}{\tau_\alpha}}$$
$$\sigma_{\beta_1} = \sqrt[2]{\frac{1}{\tau_{\beta_1}}}$$
$$\sigma_{\beta_2} = \sqrt[2]{\frac{1}{\tau_{\beta_2}}}$$

# Analysis of observed data


```r
paralleljags<-T

if(paralleljags){
    
#Source model
source("Bayesian/NmixturePoissonRagged2m.R")

#print model
writeLines(readLines("Bayesian/NmixturePoissonRagged2m.R"))

#Input Data
Dat <- c('Yobs_camera','Yobs_transect','Birds','Bird','Plant','Time','Plants','Times','resources','Nobs','cam_surveys','trans_surveys','Traitmatch')

#Inits
InitStage <- function(){
  #A blank Y matrix - all present
  initY<-array(dim=c(Birds,Plants,Times),data=1)
  list(S=initY)}

#Parameters to track
ParsStage <- c("alpha","beta1","beta2","alpha_mu","alpha_sigma","beta1_mu","beta1_sigma","beta2_mu","beta2_sigma","dtrans","dcam")

#MCMC options

ni <- 100000  # number of draws from the posterior
nt <- max(c(1,ni*.0001))  #thinning rate
nb <- ni*.9 # number to discard for burn-in
nc <- 2  # number of chains

#Jags

  Yobs_camera = mdat$Y_Camera
  Yobs_transect = mdat$Y_Transect
  Birds=max(mdat$Bird)
  Bird=mdat$Bird
  Plant=mdat$Plant
  Time=mdat$Time
  Plants=max(mdat$Plant)
  Times=max(mdat$Time)
  resources=TimeResources
  Nobs=nrow(mdat)
  cam_surveys=(mdat$Survey_Type=="Camera")*1
  trans_surveys=(mdat$Survey_Type=="Transect")*1
  Traitmatch=Traitmatch

  m<-do.call(jags.parallel,list(Dat,InitStage,ParsStage,model.file="Bayesian/NmixturePoissonRagged2m.jags",n.thin=nt, n.iter=ni,n.burnin=nb,n.chains=nc))
  
} else {
  
#Source model
source("Bayesian/NmixturePoissonRagged2m.R")

#print model
writeLines(readLines("Bayesian/NmixturePoissonRagged2m.R"))

#Input Data
Dat <- list(
  Yobs_camera = mdat$Y_Camera,
  Yobs_transect = mdat$Y_Transect,
  Birds=max(mdat$Bird),
  Bird=mdat$Bird,
  Plant=mdat$Plant,
  Time=mdat$Time,
  Plants=max(mdat$Plant),
  Times=max(mdat$Time),
  resources=TimeResources,
  Nobs=nrow(mdat),
  cam_surveys=(mdat$Y_Camera>0)*1,
  trans_surveys=(mdat$Y_Transect>0)*1,
  Traitmatch=Traitmatch)

#A blank Y matrix - all present
initY<-array(dim=c(Dat$Birds,Dat$Plants,Dat$Times),data=max(mdat$Y_Transect,na.rm=T))
initB<-as.numeric(matrix(nrow=h_species,ncol=1,data=.1))

#Inits
InitStage <- function(){list(S=initY)}

#Parameters to track
ParsStage <- c("alpha","beta1","beta2","intercept","sigma_alpha","sigma_slope1","sigma_slope2","gamma1","gamma2","dtrans","dcam")

#MCMC options

ni <- 20000  # number of draws from the posterior
nt <- max(c(1,ni*.0001))  #thinning rate
nb <- ni*.90 # number to discard for burn-in
nc <- 2  # number of chains

#Jags

m = jags(inits=InitStage,
         n.chains=nc,
         model.file="Bayesian/NmixturePoissonRagged2m.jags",
         working.directory=getwd(),
         data=Dat,
         parameters.to.save=ParsStage,
         n.thin=nt,
         n.iter=ni,
         n.burnin=nb,
         DIC=T)
}
```

```
## 
## sink("Bayesian/NmixturePoissonRagged2m.jags")
## 
## cat("
##     model {
##     #Compute true state for each pair of birds and plants
##     for (i in 1:Birds){
##     for (j in 1:Plants){
##     for (k in 1:Times){
##     
##     #Process Model
##     logit(rho[i,j,k])<-alpha[i] + beta1[i] * Traitmatch[i,j] + beta2[i] + resources[i,j,k]
##     
##     #True State
##     S[i,j,k] ~ dbern(rho[i,j,k])
##     }
##     }
##     }
##     
##     #Observation Model
##     for (x in 1:Nobs){
##     
##     #Observation Process for cameras
##     detect_cam[x]<-dcam * cam_surveys[x]
## 
##     #Observation Process for transects
##     detect_transect[x]<-dtrans * trans_surveys[x]
## 
##     Yobs_camera[x] ~ dbern(detect_cam[x] * S[Bird[x],Plant[x],Time[x]])    
##     Yobs_transect[x] ~ dbern(detect_transect[x] * S[Bird[x],Plant[x],Time[x]])    
## 
##     #Assess Model Fit - Posterior predictive check
## 
##     #Fit discrepancy statistics
##     #eval[x]<-detect[Bird[x]]*S[Bird[x],Plant[x],Camera[x]]
##     #E[x]<-pow((Yobs[x]-eval[x]),2)/(eval[x]+0.5)
##     
##     #ynew[x]~dbin(detect[Bird[x]],S[Bird[x],Plant[x],Camera[x]])
##     #E.new[x]<-pow((ynew[x]-eval[x]),2)/(eval[x]+0.5)
##     
##     }
## 
##     #Priors
##     #Observation model
##     #Detect priors, logit transformed - Following lunn 2012 p85
##     
##     #For Cameras
##     dcam ~ dunif(0,1)
## 
##     #For Transects
##     dtrans ~ dunif(0,1)
##     
##     #Detection group prior
##     #dprior_cam ~ dnorm(0,0.386)
##     #dprior_trans ~ dnorm(0,0.386)
##     
##     #Group effect detect camera
##     #tau_dcam ~ dunif(0,10)
##     #sigma_dcam<-pow(1/tau_dcam,.5)
##     
##     #Group effect detect camera
##     #tau_dtrans ~ dunif(0,100)
##     #sigma_dtrans<-pow(1/tau_dtrans,.5)
## 
##     #Process Model
##     #Species level priors
##     for (i in 1:Birds){
##       
##       #Intercept
##       #logit prior, then transform for plotting
##       alpha[i] ~ dnorm(alpha_mu,alpha_tau)
## 
##       #Traits slope 
##       beta1[i] ~ dnorm(beta1_mu,beta1_tau)    
## 
##       #Plant slope
##       beta2[i] ~ dnorm(beta2_mu,beta2_tau)    
##     }
## 
##     #Group process priors
## 
##     #Intercept 
##     alpha_mu~dnorm(0,0.386)
##     alpha_tau ~ dunif(0,100)
##     alpha_sigma<-pow(1/alpha_tau,0.5) 
##     
##     #Trait
##     beta1_mu~dnorm(0,0.386)
##     beta1_tau ~ dunif(0,100)
##     beta1_sigma<-pow(1/beta1_tau,0.5)
##     
##     #Resources
##     beta2_mu~dnorm(0,0.386)
##     beta2_tau ~ dunif(0,100)
##     beta2_sigma<-pow(1/beta2_tau,0.5)
## 
## }
##     ",fill=TRUE)
## 
## sink()
```


```r
pars<-extract_par(m)
```

###Assess Convergence


```r
###Chains
ggplot(pars[pars$par %in% c("alpha","beta1","beta2"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(par~species,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Species Level Probability")
```

![](TwoDetectSimulation_files/figure-html/unnamed-chunk-7-1.png)<!-- -->


```r
ggplot(pars[pars$par %in% c("dcam","dtrans"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + facet_grid(~par,scale="free") + theme_bw() + labs(col="Chain") + ggtitle("Detection Probability")
```

![](TwoDetectSimulation_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


```r
ggplot(pars[pars$par %in% c("alpha_mu","alpha_sigma","beta1_mu","beta1_sigma","beta2_mu","beta2_sigma"),],aes(x=Draw,y=estimate,col=as.factor(Chain))) + geom_line() + theme_bw() + labs(col="Chain") + ggtitle("Group Level Regression") + facet_wrap(~par,scales="free")
```

![](TwoDetectSimulation_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

###Posteriors


```r
###Posterior Distributions
p<-ggplot(pars[pars$par %in% c("alpha","beta1","beta2"),],aes(x=estimate)) + geom_histogram() + ggtitle("Estimate of parameters") + facet_grid(species~par,scales="free") + theme_bw() + ggtitle("Species Posteriors")

#Add true values
tr<-melt(data.frame(species=1:h_species,alpha=alpha,beta1=beta1,beta2=beta2),id.var='species')
colnames(tr)<-c("species","par","value")
psim<-p + geom_vline(data=tr,aes(xintercept=value),col='red',linetype='dashed',size=1)
#ggsave("Figures/SimulationPosteriors.jpg",dpi=300,height=8,width=8)
```


```r
p<-ggplot(pars[pars$par %in% c("beta1_mu","beta2_mu","alpha_mu","alpha_sigma","beta1_sigma","beta2_sigma","dcam","dtrans"),],aes(x=estimate)) + geom_histogram() + ggtitle("Hierarchical Posteriors") + facet_wrap(~par,scale="free",nrow=2) + theme_bw() 

#Add true values
tr<-melt(list(beta1_mu=beta1_mu,beta2_mu=beta2_mu,alpha_mu=alpha_mu,alpha_sigma=alpha_sigma,beta1_sigma=beta1_sigma,beta2_sigma=beta2_sigma,dtrans=detection_trans,dcam=detection_cam))

colnames(tr)<-c("value","par")

psim2<-p + geom_vline(data=tr,aes(xintercept=value),linetype='dashed',size=1,col="red")
#ggsave("Figures/SimulationH.jpg",dpi=300,height=4,width=10)
```

![](TwoDetectSimulation_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

###Predicted Relationship 


```r
castdf<-dcast(pars[pars$par %in% c("beta1_mu","beta2_mu","alpha_mu"),], Chain + Draw~par,value.var="estimate")

trajF<-function(alpha,beta1,beta2,x,resources){
  indat<-data.frame(alpha,beta1,beta2)
  
  #fit regression for each input estimate
  sampletraj<-list()
  
  for (y in 1:nrow(indat)){
    v=inv.logit(indat$alpha[y] + indat$beta1[y] * x + indat$beta2[y] * resources)
    
    sampletraj[[y]]<-data.frame(x=as.numeric(x),y=as.numeric(v))
  }
  
  sample_all<-rbind_all(sampletraj)
  
  #Compute CI intervals
  predy<-group_by(sample_all,x) %>% summarise(lower=quantile(y,0.025,na.rm=T),upper=quantile(y,0.975,na.rm=T),mean=mean(y,na.rm=T))
}
```

#Predicted Relationship


```r
predy<-trajF(alpha=castdf$alpha_mu,beta1=castdf$beta1_mu,x=as.numeric(traitarray),resources=as.numeric(resources),beta2=castdf$beta2_mu)

orig<-trajF(alpha=rnorm(2000,alpha_mu,alpha_sigma),beta1=rnorm(2000,beta1_mu,beta1_sigma),beta2=rnorm(2000,beta2_mu,beta2_sigma),x=as.numeric(traitarray),resources=as.numeric(resources))

#plot and compare to original data
ggplot(data=predy,aes(x=x)) + geom_point(data=mdat,aes(x=traitmatch,y=True_state),alpha=.5,size=.5)+ geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.3,fill="red")  + geom_line(aes(y=mean),size=.8,col="red",linetype="dashed") + theme_bw() + ylab("Probability of interactions") + geom_line(data=orig,aes(x=x,y=mean),col='black',size=1)+ xlab("Difference between Bill and Corolla Length") 
```

![](TwoDetectSimulation_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

The true data is plotted overtop the simulation relationship in black, and the predicted relationship in dashed red with pink CI intervals.


```r
save.image("Simulation_2M.RData")
```
