library(ggplot2)
library(dplyr)
mus<-seq(20,100,20)
thetas<-seq(20,100,20)
input<-expand.grid(mus,thetas)

out<-list()
for (x in 1:nrow(input)){
  mu<-input[x,1]
  theta<-input[x,2]
  out[[x]]<-data.frame(mu,theta,v=rnegbin(10000,mu,theta))
}
out<-dplyr::rbind_all(out)

head(out)
ggplot(out,aes(x=v,fill=as.factor(theta))) + geom_density() + facet_wrap(~mu)

out %>% group_by(mu,theta) %>% summarize(mean=mean(v),sd=sd(v))
