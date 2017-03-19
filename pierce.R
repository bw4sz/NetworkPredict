pierce<-datph %>% select(Hummingbird,Iplant_Double,ID,Month,Year,Transect_R,ele,DateP,Pierce) %>% filter(Pierce %in% c("y","Y"))


pier<-merge(pierce,flower.month,by=c("Month","Year","Transect_R"))

hpier<-pier %>% group_by(Hummingbird,R) %>% summarize(Yobs=n())
#those white-whiskered hermit observations have been invalidated, they aren't feeding.

hpier<-hpier %>% filter(!Hummingbird == "White-whiskered Hermit")
ggplot(hpier,aes(x=Hummingbird,y=Yobs,col=R)) + geom_point(size=3)  + coord_flip() + theme_bw() + labs(y="Observed piercing events",col="Resource Availability") 
ggsave("Figures/Piercing.jpeg",height=4,width=6)
