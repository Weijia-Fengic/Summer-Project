#Figure generation
rm(list=ls())
setwd("F:/EA/Project/Report")
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("maps")
library("vegan")
library("lme4")
library("minque")
library("geosphere")
library("RColorBrewer")


#load datasets :data 1 is about species-trait data ; data 2 is sample community level
data1<-read.csv("F:/EA/Project/EMP/Machado et al/EMP_Madin_gs.csv")
data2<-read.csv("F:/EA/Project/EMP/Machado et al/sample_community_summary.csv")
colourCount1<-length(unique(data1$empo_3))
getPalette1<-colorRampPalette(brewer.pal(15,"Spectral"))
#Figure 1: cell volume distribution across different habitats
#sample point
world<- map_data("world")
sample_map<-ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    color = "white", fill = "lightgray", size = 0.1) +
  geom_point(
    data = data2,
    aes(longitude_deg, latitude_deg, color = empo_3),
    alpha = 0.6,size=3.5)+
  scale_color_manual(values=getPalette1(colourCount1))+
  theme_void() +
  theme(legend.position = "none")
sample_map

#cell volume distribution
p1<-ggplot(data1,aes(x=log10(volume)))+
    geom_density(aes(fill=empo_3),lwd=1,alpha=0.8,color="black")+
    geom_histogram(aes(y=..density..),color="gray50",fill="white",alpha=0.2,position='identity')+
    xlab("Log10 cell volume(m^3)")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18))+
    theme_classic()+
      labs(fill="Habitat")+
  facet_wrap(.~empo_3,scales="free",ncol=3)+
  scale_fill_manual(values=getPalette1(colourCount1))+
  theme_bw()+
  ggtitle("Cell volume distribuion")
p1

#genome size distribution
p2<-ggplot(data1,aes(x=log10(gs_Madin)))+
  geom_density(aes(fill=empo_3),lwd=1,alpha=0.8,color="black")+
  geom_histogram(aes(y=..density..),color="gray50",fill="white",alpha=0.2,position='identity')+
  xlab("Genome size (bp)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18))+
  theme_classic()+
  labs(fill="Habitat")+
  facet_wrap(.~empo_3,scales="free",ncol=3)+
  scale_fill_manual(values=getPalette1(colourCount1))+
  theme_bw()+
  ggtitle('Genome size distribution')
p2
ggarrange(p1,p2,ncol = 1,labels = c("a","b"), common.legend = T,legend = "right")
#cell volume and genome size relationship
habitat<-split(data1,data1$empo_3)
dat<-c()
predictions<-c()
data1<-data1[!is.na(gs_Madin)]
for (i in 1:length(habitat)) {
  dat[i]<-data.frame(Genome_size=seq(from=742409.1, to=12740813.0, length=100000),Habitat=names(habitat[i]))
  predictions[i]<- predict(m2, newdata =dat[i] , type = "link", se.fit = TRUE) 
  dat[i]$pred<- predictions[i]$fit
  dat[i]$se<- predictions[i]$se.fit
  dat[i]$upperCI<- dat[i]$pred+(dat[i]$se*1.96)
  dat[i]$lowerCI<- dat[i]$pred-(dat[i]$se*1.96)
  
}  
  p4<-ggplot(dat[i], aes(x=Genome_size, y=exp(pred)))+ 
    geom_line(aes(color=Habitat))+
    geom_ribbon(aes(ymin=exp(lowerCI), ymax=exp(upperCI), fill=Habitat, alpha=0.7), show.legend = T)+ 
    geom_point(data1, mapping = aes(x=gs_Madin, y=volume, color='gray')+
    labs(y="Cell volume (m^3)", x="Genome size (bp)")+
    theme_classic()+
    scale_fill_manual(values=getPalette1(colourCount1)))+
    facet_wrap(.~empo_3,scales="free",ncol=3)
    
 
p3<-ggplot(data1,aes(x=log10(gs_Madin),y=log10(volume),fill()=empo_3))+
  geom_point(alpha=0.8,color='black')+
  geom_smooth(method='glm',color="black")+
  xlab("Log10 genome size (bp)")+
  ylab("Log10 cell volume (m^3)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18))+
  theme_classic()+
  labs(fill="Habitat")+
  facet_wrap(.~empo_3,scales="free",ncol=3)+
  scale_color_manual(values=getPalette1(colourCount1))
p3

#Shanno index
shan<-ggplot(data2,aes(x=empo_3,y=adiv_shannon,fill=empo_3))+
  geom_boxplot(color="black")+
  ylab("Shannon Index")+
  theme( axis.text.x = element_text(size=12, angle=90,vjust = 0.5,hjust = 0.6),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=15),
        axis.title.x= element_text(element_blank()),
        legend.text = element_text(size=13),
        legend.title = element_text(element_blank()))+
  theme_classic()+
  scale_color_manual(values=c('#ffff00','#8b4513','#d2b48c','#f4a460','#b8860b',
                              '#7cfc00','#006400','#00fa9a','#ffa07a','#ff6347',
                              '#ff0000','#000000','#696969','#000080','#4169e1'))
shan  



#local community in each habitat
unique(data1$empo_3)
dat1<-data1%>%filter(empo_3=="Animal surface")
dat2<-data1%>%filter(empo_3=="Animal secretion")
dat3<-data1%>%filter(empo_3=="Plant rhizosphere")
dat4<-subset(data1,empo_3=="Soil (non-saline)")
dat5<-subset(data1,empo_3=="Sediment (saline)")
dat6<-subset(data1,empo_3=="Surface (non-saline)" )
dat7<-subset(data1,empo_3=="Animal distal gut" )
dat8<-subset(data1,empo_3=="Water (non-saline)")
dat9<-subset(data1,empo_3=="Animal proximal gut")
dat10<-subset(data1,empo_3=="Water (saline)" )
dat11<-subset(data1,empo_3=="Plant corpus")
dat12<-subset(data1,empo_3=="Animal corpus")
dat13<-subset(data1,empo_3=="Plant surface" )
dat14<-subset(data1,empo_3=="Sediment (non-saline)")
dat15<-subset(data1,empo_3=="Surface (saline)")

#community composition 
com1<-dat1 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com1<-pivot_wider(com1,names_from = "Species",values_from = "n")
com1[is.na(com1)]<-0
rownames(com1)<-com1$sample
com1<-com1[,-1]
com1_abun<-decostand(com1,method="total")
com1_dis<-vegdist(com1_abun,method = "bray")
v1<-data2%>%filter(empo_3=="Animal surface")%>%select(V_mean,V_var,V_median,V_max,V_min)
v1_dist<-dist(v1,method='euclidean')
mantel(com1_dis,v1_dist,method='spearman')

com2<-dat2 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com2<-pivot_wider(com2,names_from = "Species",values_from = "n")
com2[is.na(com2)]<-0
rownames(com2)<-com2$sample
com2<-com2[,-1]
com2_abun<-decostand(com2,method="total")
com2_dis<-vegdist(com2_abun,method = "bray")
v2<-data2%>%filter(empo_3=="Animal secretion")%>%select(V_mean,V_var,V_median,V_max,V_min)
v2_dist<-dist(v2,method='euclidean')
mantel(com2_dis,v2_dist,method='spearman')
mantel(dist.geo2,v2_dist,method = 'spearman')
mantel(disT2,v2_dist,method = 'spearman')
mantel(disph2,v2_dist,method = 'spearman')

com3<-dat3 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com3<-pivot_wider(com3,names_from = "Species",values_from = "n")
com3[is.na(com3)]<-0
rownames(com3)<-com3$sample
com3<-com3[,-1]
com3_abun<-decostand(com3,method="total")
com3_dis<-vegdist(com3_abun,method = "bray")
v3<-data2%>%filter(empo_3=="Plant rhizosphere")%>%select(V_mean,V_var,V_median,V_max,V_min)
v3_dist<-dist(v3,method='euclidean')
mantel(com3_dis,v3_dist,method='spearman')
mantel(dist.geo3,v3_dist,method = 'spearman')

com4<-dat4 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com4<-pivot_wider(com4,names_from = "Species",values_from = "n")
com4[is.na(com4)]<-0
rownames(com4)<-com4$sample
com4<-com4[,-1]
com4_abun<-decostand(com4,method="total")
com4_dis<-vegdist(com4_abun,method = "bray")
v4<-data2%>%filter(empo_3=="Soil (non-saline)")%>%select(V_mean,V_var,V_median,V_max,V_min)
v4_dist<-dist(v4,method='euclidean')
mantel(com4_dis,v4_dist,method='spearman')
mantel(dist.geo4,v4_dist,method = 'spearman')


com5<-dat5 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com5<-pivot_wider(com5,names_from = "Species",values_from = "n")
com5[is.na(com5)]<-0
rownames(com5)<-com5$sample
com5<-com5[,-1]
com5_abun<-decostand(com5,method="total")
com5_dis<-vegdist(com5_abun,method = "bray")
v5<-data2%>%filter(empo_3=="Sediment (saline)")%>%select(V_mean,V_var,V_median,V_max,V_min)
v5_dist<-dist(v5,method='euclidean')
mantel(com5_dis,v5_dist,method='spearman')
mantel(dist.geo5,v5_dist,method = 'spearman')

com6<-dat6 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com6<-pivot_wider(com6,names_from = "Species",values_from = "n")
com6[is.na(com6)]<-0
rownames(com6)<-com6$sample
com6<-com6[,-1]
com6_abun<-decostand(com6,method="total")
com6_dis<-vegdist(com6_abun,method = "bray")
v6<-data2%>%filter(empo_3=="Surface (non-saline)")%>%select(V_mean,V_var,V_median,V_max,V_min)
v6_dist<-dist(v6,method='euclidean')
mantel(com6_dis,v6_dist,method='spearman')
mantel(dist.geo6,v6_dist,method = 'spearman')

com7<-dat7 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com7<-pivot_wider(com7,names_from = "Species",values_from = "n")
com7[is.na(com7)]<-0
rownames(com7)<-com7$sample
com7<-com7[,-1]
com7_abun<-decostand(com7,method="total")
com7_dis<-vegdist(com7_abun,method = "bray")
v7<-data2%>%filter(empo_3=="Animal distal gut")%>%select(V_mean,V_var,V_median,V_max,V_min)
v7_dist<-dist(v7,method='euclidean')
mantel(com7_dis,v7_dist,method='spearman')
mantel(dist.geo7,v7_dist,method = 'spearman')

com8<-dat8 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com8<-pivot_wider(com8,names_from = "Species",values_from = "n")
com8[is.na(com8)]<-0
rownames(com8)<-com8$sample
com8<-com8[,-1]
com8_abun<-decostand(com8,method="total")
com8_dis<-vegdist(com8_abun,method = "bray")
v8<-data2%>%filter(empo_3=="Water (non-saline)")%>%select(V_mean,V_var,V_median,V_max,V_min)
v8_dist<-dist(v8,method='euclidean')
mantel(com8_dis,v8_dist,method='spearman')
mantel(dist.geo8,v8_dist,method = 'spearman')

com9<-dat9 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com9<-pivot_wider(com9,names_from = "Species",values_from = "n")
com9[is.na(com9)]<-0

com10<-dat10 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com10<-pivot_wider(com10,names_from = "Species",values_from = "n")
com10[is.na(com10)]<-0

com11<-dat11 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com11<-pivot_wider(com11,names_from = "Species",values_from = "n")
com11[is.na(com11)]<-0

com12<-dat12 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com12<-pivot_wider(com12,names_from = "Species",values_from = "n")
com12[is.na(com12)]<-0

com13<-dat13 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com13<-pivot_wider(com13,names_from = "Species",values_from = "n")
com13[is.na(com13)]<-0

com14<-dat14 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com14<-pivot_wider(com14,names_from = "Species",values_from = "n")
com14[is.na(com14)]<-0

com15<-dat15 %>%
  group_by(sample, Species) %>%
  summarise(n =n(), .groups = "drop") %>%
  filter(n >= 1L) 
com15<-pivot_wider(com15,names_from = "Species",values_from = "n")
com15[is.na(com15)]<-0

#environmental distance matrix
##pH matrix
ph<-data2%>%filter(ph!='NA')
Ph1<-data2%>%filter(empo_3=="Animal surface")%>%select(ph)
disph1<-dist(Ph1,method="euclidean")
Ph2<-data2%>%filter(empo_3=="Animal secretion")%>%select(ph)
disph2<-dist(Ph2,method="euclidean")
Ph3<-data2%>%filter(empo_3=="Plant rhizosphere")%>%select(ph)
disph3<-dist(Ph3,method="euclidean")

Ph4<-data2%>%filter(empo_3=="Soil (non-saline)")%>%select(ph)
disph4<-dist(Ph4,method="euclidean")
Ph5<-data2%>%filter(empo_3=="Sediment (saline)" )%>%select(ph)
disph5<-dist(Ph5,method="euclidean")
Ph6<-data2%>%filter(empo_3=="Surface (non-saline)")%>%select(ph)
disph6<-dist(Ph6,method="euclidean")

Ph7<-data2%>%filter(empo_3=="Animal distal gut")%>%select(ph)
disph7<-dist(Ph7,method="euclidean")

Ph8<-ph%>%filter(empo_3=="Water (non-saline)")%>%select(ph)
disph8<-dist(Ph8,method="euclidean")
v8<-ph%>%filter(empo_3=="Water (non-saline)")%>%select(V_mean,V_var,V_median,V_max,V_min)
v8_dist<-dist(v8,method='euclidean')
mantel(disph8,v8_dist,method='pearman')


Ph9<-data2%>%filter(empo_3=="Animal proximal gut")%>%select(ph)
disph9<-dist(Ph9,method="euclidean")

Ph10<-data2%>%filter(empo_3=="Water (saline)" )%>%select(ph)
disph10<-dist(Ph10,method="euclidean")
Ph11<-data2%>%filter(empo_3=="Plant corpus")%>%select(ph)
disph11<-dist(Ph11,method="euclidean")
Ph12<-data2%>%filter(empo_3=="Animal corpus" )%>%select(ph)
disph12<-dist(Ph12,method="euclidean")

Ph13<-data2%>%filter(empo_3=="Plant surface" )%>%select(ph)
disph13<-dist(Ph13,method="euclidean")
Ph14<-data2%>%filter(empo_3=="Sediment (non-saline)")%>%select(ph)
disph14<-dist(Ph14,method="euclidean")
Ph15<-data2%>%filter(empo_3=="Surface (saline)")%>%select(ph)
disph15<-dist(Ph15,method="euclidean")

#temp distance
temp<-data2%>%filter(temperature_deg_c!='NA')
temp1<-data2%>%filter(empo_3=="Animal surface")%>%select(temperature_deg_c)
disT1<-dist(temp1,method="euclidean")
temp2<-data2%>%filter(empo_3=="Animal secretion")%>%select(temperature_deg_c)
disT2<-dist(temp2,method="euclidean")
temp3<-data2%>%filter(empo_3=="Plant rhizosphere")%>%select(temperature_deg_c)
disT3<-dist(temp3,method="euclidean")

temp4<-data2%>%filter(empo_3=="Soil (non-saline)")%>%select(temperature_deg_c)
disT4<-dist(temp4,method="euclidean")
temp5<-data2%>%filter(empo_3=="Sediment (saline)" )%>%select(temperature_deg_c)
disT5<-dist(temp5,method="euclidean")
temp6<-data2%>%filter(empo_3=="Surface (non-saline)")%>%select(temperature_deg_c)
disT6<-dist(temp6,method="euclidean")

temp7<-data2%>%filter(empo_3=="Animal distal gut")%>%select(temperature_deg_c)
disT7<-dist(temp7,method="euclidean")

temp8<-temp%>%filter(empo_3=="Water (non-saline)")%>%select(temperature_deg_c)
disT8<-dist(temp8,method="euclidean")
v8<-temp%>%filter(empo_3=="Water (non-saline)")%>%select(V_mean,V_var,V_median,V_max,V_min)
v8_dist<-dist(v8,method='euclidean')
mantel(disT8,v8_dist,method='spearman')

temp9<-data2%>%filter(empo_3=="Animal proximal gut")%>%select(temperature_deg_c)
disT9<-dist(temp9,method="euclidean")

temp10<-data2%>%filter(empo_3=="Water (saline)" )%>%select(temperature_deg_c)
disT10<-dist(temp10,method="euclidean")
temp11<-data2%>%filter(empo_3=="Plant corpus")%>%select(temperature_deg_c)
disT11<-dist(temp11,method="euclidean")
temp12<-data2%>%filter(empo_3=="Animal corpus" )%>%select(temperature_deg_c)
disT12<-dist(temp12,method="euclidean")

temp13<-temp%>%filter(empo_3=="Plant surface" )%>%select(temperature_deg_c)
disT13<-dist(temp13,method="euclidean")
v13<-temp%>%filter(empo_3=="Plant surface")%>%select(V_mean,V_var,V_median,V_max,V_min)
v13_dist<-dist(v13,method='euclidean')
mantel(disT13,v13_dist,method='spearman')

temp14<-data2%>%filter(empo_3=="Sediment (non-saline)")%>%select(temperature_deg_c)
disT14<-dist(temp14,method="euclidean")
temp15<-data2%>%filter(empo_3=="Surface (saline)")%>%select(temperature_deg_c)
disT15<-dist(temp15,method="euclidean")


###geographic distance
geo1<-data2%>%filter(empo_3=="Animal surface")%>%select(longitude_deg,latitude_deg)
d.geo1<-distm(geo1,fun=distHaversine)
dist.geo1<-as.dist(d.geo1)

geo2<-data2%>%filter(empo_3=="Animal secretion" )%>%select(longitude_deg,latitude_deg)
d.geo2<-distm(geo2,fun=distHaversine)
dist.geo2<-as.dist(d.geo2)

geo3<-data2%>%filter(empo_3=="Plant rhizosphere")%>%select(longitude_deg,latitude_deg)
d.geo3<-distm(geo3,fun=distHaversine)
dist.geo3<-as.dist(d.geo3)

geo4<-data2%>%filter(empo_3=="Soil (non-saline)")%>%select(longitude_deg,latitude_deg)
d.geo4<-distm(geo4,fun=distHaversine)
dist.geo4<-as.dist(d.geo4)

geo5<-data2%>%filter(empo_3=="Sediment (saline)")%>%select(longitude_deg,latitude_deg)
d.geo5<-distm(geo5,fun=distHaversine)
dist.geo5<-as.dist(d.geo5)

geo6<-data2%>%filter(empo_3=="Surface (non-saline)" )%>%select(longitude_deg,latitude_deg)
d.geo6<-distm(geo6,fun=distHaversine)
dist.geo6<-as.dist(d.geo6)

geo7<-data2%>%filter(empo_3=="Animal distal gut")%>%select(longitude_deg,latitude_deg)
d.geo7<-distm(geo7,fun=distHaversine)
dist.geo7<-as.dist(d.geo7)

geo8<-data2%>%filter(empo_3=="Water (non-saline)")%>%select(longitude_deg,latitude_deg)
d.geo8<-distm(geo8,fun=distHaversine)
dist.geo8<-as.dist(d.geo8)

geo9<-data2%>%filter(empo_3=="Animal proximal gut")%>%select(longitude_deg,latitude_deg)
d.geo9<-distm(geo9,fun=distHaversine)
dist.geo9<-as.dist(d.geo9)

geo10<-data2%>%filter(empo_3=="Water (saline)")%>%select(longitude_deg,latitude_deg)
d.geo10<-distm(geo10,fun=distHaversine)
dist.geo10<-as.dist(d.geo10)

geo11<-data2%>%filter(empo_3=="Plant corpus" )%>%select(longitude_deg,latitude_deg)
d.geo11<-distm(geo11,fun=distHaversine)
dist.geo11<-as.dist(d.geo11)

geo12<-data2%>%filter(empo_3=="Animal corpus" )%>%select(longitude_deg,latitude_deg)
d.geo12<-distm(geo12,fun=distHaversine)
dist.geo12<-as.dist(d.geo12)

geo13<-data2%>%filter(empo_3=="Plant surface")%>%select(longitude_deg,latitude_deg)
d.geo13<-distm(geo13,fun=distHaversine)
dist.geo13<-as.dist(d.geo13)

geo14<-data2%>%filter(empo_3=="Sediment (non-saline)")%>%select(longitude_deg,latitude_deg)
d.geo14<-distm(geo14,fun=distHaversine)
dist.geo14<-as.dist(d.geo14)

geo15<-data2%>%filter(empo_3=="Surface (saline)")%>%select(longitude_deg,latitude_deg)
d.geo15<-distm(geo15,fun=distHaversine)
dist.geo15<-as.dist(d.geo15)

#Mantel test
##mean volume VS composition 
mmph1<-mantel(disvmean1,comdist1,method='spearman',permutations = 999,na.rm=TRUE)
mmph2<-mantel(disvmean2,comdist2,method='spearman',permutations = 9999,na.rm=TRUE)
mmph3<-mantel(disvmean3,comdist3,method='spearman',permutations = 9999,na.rm=TRUE)
mmph4<-mantel(disvmean4,comdist4,method='spearman',permutations = 9999,na.rm=TRUE)
mmph6<-mantel(disvmean5,comdist5,method='spearman',permutations = 9999,na.rm=TRUE)
mmph7<-mantel(disvmean7,comdist7,method='spearman',permutations = 9999,na.rm=TRUE)
mmph8<-mantel(disvmean8,comdist8,method='spearman',permutations = 9999,na.rm=TRUE)
mmph9<-mantel(disvmean9,comdist9,method='spearman',permutations = 9999,na.rm=TRUE)
mmph10<-mantel(disvmean10,comdist10,method='spearman',permutations = 9999,na.rm=TRUE)
mmph11<-mantel(disvmean11,comdist11,method='spearman',permutations = 9999,na.rm=TRUE)
mmph12<-mantel(disvmean12,comdist12,method='spearman',permutations = 9999,na.rm=TRUE)
mmph13<-mantel(disvmean13,comdist13,method='spearman',permutations = 9999,na.rm=TRUE)
mmph14<-mantel(disvmean14,comdist14,method='spearman',permutations = 9999,na.rm=TRUE)
mmph15<-mantel(disvmean15,comdist15,method='spearman',permutations = 9999,na.rm=TRUE)

#mean volume VS geographical distance
mmgeo1<-mantel(disvmean1,dist.geo1,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo2<-mantel(disvmean2,dist.geo2,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo3<-mantel(disvmean3,dist.geo3,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo4<-mantel(disvmean4,dist.geo4,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo5<-mantel(disvmean5,dist.geo5,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo6<-mantel(disvmean6,dist.geo6,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo7<-mantel(disvmean7,dist.geo7,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo8<-mantel(disvmean8,dist.geo8,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo9<-mantel(disvmean9,dist.geo9,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo10<-mantel(disvmean10,dist.geo10,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo11<-mantel(disvmean11,dist.geo11,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmgeo12<-mantel(disvmean12,dist.geo12,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo13<-mantel(disvmean13,dist.geo13,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo14<-mantel(disvmean14,dist.geo14,method = 'spearman',permutations = 999,na.rm = TRUE)
mmgeo15<-mantel(disvmean15,dist.geo15,method = 'spearman',permutations = 999,na.rm = TRUE)

#mean volume VS temp
mmT1<-mantel(disvmean1,disT1,method = 'spearman',permutations = 999,na.rm = TRUE)
mmT2<-mantel(disvmean2,distT2,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT3<-mantel(disvmean3,distT3,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT4<-mantel(disvmean4,distT4,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT5<-mantel(disvmean5,distT5,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT6<-mantel(disvmean6,distT6,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT7<-mantel(disvmean7,distT7,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT8<-mantel(disvmean8,distT8,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT9<-mantel(disvmean9,distT9,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT10<-mantel(disvmean10,distT10,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT11<-mantel(disvmean11,distT11,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT12<-mantel(disvmean12,distT12,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT13<-mantel(disvmean13,distT13,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT14<-mantel(disvmean14,distT14,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmT15<-mantel(disvmean15,distT15,method = 'spearman',permutations = 9999,na.rm = TRUE)

#mean volume VS pH
mmph1<-mantel(disvmean1,disph1,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph2<-mantel(disvmean2,disph2,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph3<-mantel(disvmean3,disph3,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph4<-mantel(disvmean4,disph4,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph5<-mantel(disvmean5,disph5,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph6<-mantel(disvmean6,disph6,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph7<-mantel(disvmean7,disph7,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph8<-mantel(disvmean8,disph8,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph9<-mantel(disvmean9,disph9,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph10<-mantel(disvmean10,disph10,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph11<-mantel(disvmean11,disph11,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph12<-mantel(disvmean12,disph12,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph13<-mantel(disvmean13,disph13,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph14<-mantel(disvmean14,disph14,method = 'spearman',permutations = 9999,na.rm = TRUE)
mmph15<-mantel(disvmean15,disph15,method = 'spearman',permutations = 9999,na.rm = TRUE)

##visualization 
r<-read.csv("F:/EA/Project/Report/comparing R^2 in Madin and Machado.csv")
ggplot(r,aes(x=Habitat,y=R2,color=Habitat,shape=Data.source))+
  geom_point(size=4,alpha=0.9)+
  geom_point(size=4,alpha=0.9)+
  theme(axis.text = element_text(size=13),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size = 13, angle = 90, vjust = 0.5, hjust = 1))+
  theme_classic()+
  labs(fill="Habitat",
       x="",
       y="R^2")+
  scale_color_manual(values=getPalette1(colourCount1))
  
  


