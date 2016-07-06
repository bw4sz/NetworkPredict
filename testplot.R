predy<-trajF(alpha=-1.43,beta1=-0.056,trait=indat$Traitmatch,resources=scale(indat$All_Flowers),beta2=-0.11,beta3=0.001)
predy2<-trajF(alpha=-1.43,beta1=-0.056,trait=indat$Traitmatch,resources=scale(indat$All_Flowers),beta2=0,beta3=0.001)

ggplot(data=predy,aes(x=trait)) + geom_ribbon(aes(ymin=lower,ymax=upper),alpha=0.4,fill="red")  +  theme_bw() + ylab("Interactions") + xlab("Difference between Bill and Corolla Length")  + geom_line(aes(y=mean)) + geom_line(data=predy2,aes(x=trait,y=mean),col="blue")

#Trajectories from posterior
predH<-trajF(alpha=3.4066,beta1=-0.15,trait=indat[indat$BAll_Flowers==1,"Traitmatch"],resources=indat[indat$ BAll_Flowers==1,"scaledR"],beta2=-0.4,beta3=0.04)

predL<-trajF(alpha=3.4066,beta1=-0.15,trait=indat[indat$BAll_Flowers==0,"Traitmatch"],resources=indat[indat$ BAll_Flowers==0,"scaledR"],beta2=-0.4,beta3=0.1)

predhl<-melt(list(High=predH,Low=predL),id.vars=colnames(predH))

colnames(predhl)[5]<-"BFlowerL"

indat$BFlowerL<-factor(as.character(indat$BAll_Flowers))
levels(indat$BFlowerL)<-c("Low","High")

ggplot(data=predhl,aes(x=trait)) + geom_ribbon(aes(ymin=lower,ymax=upper,fill=BFlowerL),alpha=0.2)  + geom_line(aes(y=mean,col=BFlowerL),size=.8) + theme_bw() + ylab("Interactions") + xlab("Difference between Bill and Corolla Length") + geom_point(data=mindat,aes(x=Traitmatch,y=value))+ labs(fill="Resource Availability",col="Resource Availability") 
