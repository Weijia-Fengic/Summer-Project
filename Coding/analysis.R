#Statistical analysis and model fitting 
rm(list=ls())
setwd("F:/EA/Project/Report")
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("vegan")
library("lme4")
library("minque")
library("geosphere")
library("coin")
library("moments")
library("car")
library("AICcmodavg")
library("MuMIn")

#load datasets :data 1 is about species-trait data ; data 2 is sample community level
#replace 'NA' in gs_Madin from gs_Machado
data1<-read.csv("F:/EA/Project/EMP/Machado et al/EMP_Madin_gs.csv")
data1<-data1%>%mutate(gs_Madin=coalesce(gs_Madin,gs_Machado))
data2<-read.csv("F:/EA/Project/Report/sample_community_summary.csv")

#descriptive statistics
#overall
log10(mean(data1$volume))
log10(var(data1$volume))
log10(range(data1$volume))
log10(median(data1$volume))

log10(mean(data1$gs_Madin))
log10(var(data1$gs_Madin))
log10(range(data1$gs_Madin))
log10(median(data1$gs_Madin))
v<-ggplot(data1,aes(x=empo_3,y=log10(volume),fill=empo_3))+
  geom_boxplot(color='black')+
  labs(x="",y="Log10 cell volume",fill="Habitat")+
  theme_bw()+
  theme_classic()+
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size=15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size=15))+
  scale_fill_manual(values=getPalette1(colourCount1))

gs<-ggplot(data1,aes(x=empo_3,y=log10(gs_Madin),fill=empo_3))+
  geom_boxplot(color='black')+
  labs(x="",y="Log10 genome size",fill="Habitat")+
  theme_bw()+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size=15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size=15))+
  scale_fill_manual(values=getPalette1(colourCount1))

ggarrange(v,gs,labels = c("a","b"),ncol=1,common.legend = T,legend = "right")

#in 15 habitats 
general<-data1%>%group_by(empo_3)%>%
  summarise(mean=log10(mean(volume)),
            median=log10(median(volume)),
            var=log10(var(volume)),
            max=log10(max(volume)),
            min=log10(min(volume)))

general1<-data1%>%group_by(empo_3)%>%
  summarise(mean=log10(mean(gs_Madin)),
            median=log10(median(gs_Madin)),
            var=log10(var(gs_Madin)),
            max=log10(max(gs_Madin)),
            min=log10(min(gs_Madin)))

data1$x<-log10(data1$gs_Madin)
data1$y<-log10(data1$volume)
data1$z<-log10(data1$gs_Machado)
Coef2<-data1 %>% group_by(empo_3) %>% 
  summarise(
    r.sqr = summary(lm(y~z))$r.squared,
    Intercept=summary(lm(y~z))$coefficients[[1]],
    slop = summary(lm(y~z))$coefficients[[2]],
    p_value=summary(lm(y~z))$coefficients[[3]])
newdata1<-merge(Coef1,data,by="empo_3")
newdata2<-newdata1%>%
  group_by(empo_3)%>%
  summarise(R2=sum(r.sqr),
            intercept=sum(Intercept),
            Slop=sum(slop),
            Pr=sum(p_value))

r<-cbind.data.frame(Coef1$empo_3,Coef1$r.sqr,Coef2$r.sqr)
names(r)[1:3]<-c("Habitat","Madin","Machado")
write.csv(r,"comparing R^2 in Madin and Machado.csv")
r_compare<-ggplot(r,aes(x=Habitat,color=Habitat,fill=Habitat))+
  geom_point(aes(y=Madin),size=4,shape=2,alpha=0.8)+
  geom_point(aes(y=Machado),size=4,shape=1,alpha=0.8)+
  labs(x="",y="R^2")+
theme_classic()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size=15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size=15))
  

#species with largest genome size 
sp_name<-c()
group<-c()
for (i in 1:nrow(data1)) {
  if(data1[i,'gs_Madin']==585297){
    sp_name<-c(sp_name,data1[i,'Species'])
    group<-c(group,data1[i,'empo_3'])
  }
}
unique(sp_name)
unique(group)


#Normal distribution: Q-Q plots
habitat<-split(data1,data1$empo_3)
plot_list<-list()
pdf("Q-Q plots.pdf")
for (i in 1:length(habitat)) {
  sub_data1<-habitat[[i]]
  n1<-names(habitat[i])
  p1<-ggqqplot(log10(sub_data1$volume))+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20))+
    theme_classic()+
    labs(title=names(habitat[i]))
  plot_list[[i]]<-p1 
  print(plot_list[[i]])
  
}
dev.off()

#Normal distribution: skewness and kurtosis test
p.value<-c()
env<-c()
for (i in 1:length(habitat)) {
 dat<-habitat[[i]]
  n<-jarque.test(log10(dat$volume))
  env<-c(env,names(habitat[i]))
  p.value<-c(p.value,n$p.value)
}
nor<-cbind.data.frame(env,p.value)

  
#Explore whether there is difference in cell volume distribution across different habitats
#Test for variance
leveneTest(log10(volume)~empo_3,data=data1)
leveneTest(log10(gs_Madin)~empo_3,data=data1)
#Test for mean
test<-aov(volume~empo_3,data=data1)
summary(test)
TukeyHSD(test)

summary(aov(gs_Madin~empo_3,data=data1))
m1<-glm(volume~gs_Madin*empo_3,data=data1)
summary(m1)
anova(m1,test="Chisq")

#Phylum Abundance table
com<-data1%>%
  group_by(sample,Species)%>%
  summarise(n=n())%>%
  filter(n>=1L)
com_abun<-pivot_wider(com,names_from = "Species",values_from = "n")
com_abun[is.na(com_abun)]<-0
head(rownames(com_abun))
head(colnames(com_abun))
rownames(com_abun)<-com_abun$sample
com_abun<-com_abun[,-1]
com_abun<-decostand(com_abun,method="total")
write.csv(com_abun,"community species relative abundance.csv ")
#calculate relative abundance in each habitat
com<-data1%>%group_by(empo_3,phylum)%>%
  summarise(n=n())%>%
  filter(n>=1L)
com_abun<-pivot_wider(com,names_from = "phylum",values_from = "n")
com_abun[is.na(com_abun)]<-0
head(rownames(com_abun))
head(colnames(com_abun))
rownames(com_abun)<-com_abun$empo_3
com_abun<-com_abun[,-1]
com_abun<-decostand(com_abun,method="total")
com_abun$empo_3<-unique(data1$empo_3)
com_phy<-com_abun%>%
 pivot_longer(names_to = "phylum",values_to = "Abundance",
              cols = !empo_3)
com_phy$abundance<-(com_phy$Abundance)*100
phy_habitat<-ggplot(com_phy,aes(x=empo_3,y=abundance,fill=phylum))+
  geom_bar(stat = "identity")+
  labs(x="",y="Relative Abundance (%)")+
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size=15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size=15))
phy_habitat
pv<-aov(volume~phylum,data=data1)
summary(aov(gs_Madin~phylum,data=data1))

#Community composition 
species<-read.csv("F:/EA/Project/Report/community species relative abundance.csv")
head(row.names(species))
species<-species[,-1]
row.names(species)<-data2$sample
#PCoa 
distance<-vegdist(species,method = "bray")
com_distance<-as.matrix(distance)
write.csv(com_distance,"sample composition dismiliraty.csv ")
pcoa<-cmdscale(distance,k=(nrow(species)-1),eig=TRUE)

#Pcoa visualization 
plot_data<-data.frame({pcoa$points})[1:3]
plot_data$sample<-rownames(plot_data)
names(plot_data)[1:3]<-c("PCoA1","PCoA2","PCoA3")
eig<-pcoa$eig
plot_data$Habitat<-data2$empo_3
summary(aov(PCoA1~Habitat,data=plot_data))
summary(aov(PCoA2~Habitat,data=plot_data))
write.csv(plot_data,"PCoA_plot.csv")
com_com<-ggplot(data=plot_data,aes(PCoA1,PCoA2,color=Habitat))+
  geom_point(alpha=0.7,size=2)+
  stat_chull(fill=NA)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  theme_classic()+
  scale_color_manual(values=c('#ffff00','#8b4513','#d2b48c','#f4a460','#b8860b',
                              '#7cfc00','#006400','#00fa9a','#ffa07a','#ff6347',
                              '#ff0000','#000000','#696969','#000080','#4169e1'))+
  theme(axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size=15),
    axis.title.x= element_text(size=15),
    legend.text = element_text(size = 13),
    legend.title = element_text(size=15))

com_com1<-ggplot(data=plot_data,aes(PCoA1,PCoA3,color=Habitat))+
  geom_point(alpha=0.7,size=2)+
  stat_chull(fill=NA)+
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""))+
  theme_classic()+
  scale_color_manual(values=c('#ffff00','#8b4513','#d2b48c','#f4a460','#b8860b',
                              '#7cfc00','#006400','#00fa9a','#ffa07a','#ff6347',
                              '#ff0000','#000000','#696969','#000080','#4169e1'))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size=15),
        axis.title.x= element_text(size=15),
        legend.text = element_text(size = 13),
        legend.title = element_text(size=15))

#model fitting
plot_data<-read.csv("F:/EA/Project/Report/PCoA_plot.csv")
data2$PCoA1<-plot_data$PCoA1
data2$PCoA2<-plot_data$PCoA2
data2$PcoA3<-plot_data$PCoA3
data2$V_mean<-scale(data2$V_mean)
data2$V_var<-scale(data2$V_var)
data2$GS_mean.y<-scale(data2$GS_mean.y)
data2$GS_var.y<-scale(data2$GS_var.y)
data2$lat<-scale(data2$latitude_deg)
data2$long<-scale(data2$longitude_deg)
data2$ph<-scale(data2$ph)
data2$temp<-scale(data2$temperature_deg_c)
data2$DO<-scale(data2$oxygen_mg_per_l)
data2$depth<-scale(data2$depth_m)
data2$altitude<-scale(data2$altitude_m)
data2$elevation<-scale(data2$elevation_m)
data2$P<-scale(data2$phosphate_umol_per_l)
data2$NH<-scale(data2$ammonium_umol_per_l)
data2$NO<-scale(data2$nitrate_umol_per_l)
data2$S<-scale(data2$sulfate_umol_per_l)
m1<-lmer(V_mean~GS_mean.y+PCoA1+PCoA2+lat+long+ph+temp+depth+altitude+elevation+(1|empo_3),data=data2)
summary(m1)
subdata2<-data2[,!is.na('V_var')]
m2<-lmer(V_var~GS_var.y+PCoA1+PCoA2+lat+long+ph+temp+depth+altitude+elevation+(1|empo_3),data=subdata2)
summary(m2)
#Pick the most important varialbles
###generate global modal and find fitness model
options(na.action = "na.fail") #Required for dredge to run
model_dredge2<-dredge(m2,beta=F,evaluate=T,rank=AICc)
options(na.action = 'na.omit') #set back to default 
head(model_dredge2)
nrow(model_dredge) #the number of model 
top_model2<-get.models(model_dredge2,subset=1)[[1]]#find the top model and interpret it
top_model2
get.models(model_dredge3,subset=delta<=2)
get.models(model_dredge3,subset=1)
###evidence to support top model
summary(model.avg(model_dredge,subset=delta<2))
###Summarize the top model 
summary(top_model2) #interpret this top model 
vif(top_model2)#test the robustness of the model
vif(model_dredge,subset=delta<=2)



# trait dissimilarity explained among habitats
volume_dist<-dist(data2[,2:6],method='euclidean')
gs_dist<-dist(data2[,24:28],method='euclidean')
volume<-as.vector(volume_dist)
gs<-as.factor(gs_dist)
anosim(volume_dist~ empo_3, data = data2)

#environmental factors
#dissolve oxygen
do<-filter(data2,oxygen_mg_per_l!='NA')
unique(do$empo_3)
v_do1<-do%>%filter(empo_3=='Water (non-saline)')%>%select(V_mean,V_var,V_median,V_max,V_min)
v1<-dist(v_do1,method='euclidean')
d1<-do%>%filter(empo_3=='Water (non-saline)')
d1.dist<-dist(d1,method='euclidean')
mantel(v1,d1.dist,method = 'spearman')

#Phosphorus
phos<-filter(data2,phosphate_umol_per_l!='NA')
unique(phos$empo_3)
p1<-filter(phos,empo_3=='Water (saline)')
p1.dist<-dist(p1,method='euclidean')
p2<-filter(phos,empo_3=='Water (non-saline)')
p2.dist<-dist(p2,method='euclidean')
p3<-filter(phos,empo_3=='Sediment (saline)')
p3.dist<-dist(p3,method='euclidean')
p4<-filter(phos,empo_3=='Surface (saline)')
p4.dist<-dist(p4,method='euclidean')
vp1<-phos%>%filter(empo_3=='Water (saline)')%>%select(V_mean,V_var,V_median,V_max,V_min)
v1<-dist(vp1,method='euclidean')
vp2<-phos%>%filter(empo_3=='Water (non-saline)')%>%select(V_mean,V_var,V_median,V_max,V_min)
v2<-dist(vp2,method='euclidean')
vp3<-phos%>%filter(empo_3=='Sediment (saline)')%>%select(V_mean,V_var,V_median,V_max,V_min)
v3<-dist(vp3,method='euclidean')
vp4<-phos%>%filter(empo_3=='Surface (saline)')%>%select(V_mean,V_var,V_median,V_max,V_min)
v4<-dist(vp4,method='euclidean')
mantel(v1,p1.dist,method = 'spearman')
mantel(v2,p2.dist,method = 'spearman')
mantel(v3,p3.dist,method = 'spearman')
mantel(v4,p4.dist,method = 'spearman')

#Ammonium
NH<-filter(data2,ammonium_umol_per_l!='NA')
unique(NH$empo_3)
nh1<-filter(NH,empo_3=='Water (saline)')
nh1.dist<-dist(nh1$ammonium_umol_per_l,method='euclidean')
vnh1<-dist(nh1[,3:7],method='euclidean')
nh2<-filter(NH,empo_3=='Soil (non-saline)')
nh2.dist<-dist(nh2$ammonium_umol_per_l,method='euclidean')
vnh2<-dist(nh2[,3:7],method='euclidean')
nh3<-filter(NH,empo_3=='Water (non-saline)')
nh3.dist<-dist(nh3$ammonium_umol_per_l,method='euclidean')
vnh3<-dist(nh3[,3:7],method='euclidean')
nh4<-filter(NH,empo_3=='Sediment (saline)')
nh4.dist<-dist(nh4$ammonium_umol_per_l,method='euclidean')
vnh4<-dist(nh4[,3:7],method='euclidean')
mantel(vnh1,nh1.dist,method = 'spearman')
mantel(vnh2,nh2.dist,method = 'spearman')
mantel(vnh3,nh3.dist,method = 'spearman')
mantel(vnh4,nh4.dist,method = 'spearman')

#Nitrate
NO<-filter(data2,nitrate_umol_per_l!='NA')
unique(NO$empo_3)
n1<-filter(NO,empo_3=='Water (saline)')
n2<-filter(NO,empo_3=='Soil(non-saline)')
n3<-filter(NO,empo_3=='Water (non-saline)')
n4<-filter(NO,empo_3=='Sediment (saline)')
vn1<-dist(n1[,3:7],method='euclidean')
vn3<-dist(n3[,3:7],method='euclidean')
vn4<-dist(n4[,3:7],method='euclidean')
n1.dist<-dist(n1$nitrate_umol_per_l,method = 'euclidean')
n3.dist<-dist(n3$nitrate_umol_per_l,method = 'euclidean')
n4.dist<-dist(n4$nitrate_umol_per_l,method = 'euclidean')
mantel(vn1,n1.dist,method='spearman')
mantel(vn3,n3.dist,method='spearman')
mantel(vn4,n4.dist,method='spearman')

#Sulfate
S<-filter(data2,sulfate_umol_per_l!='NA')
unique(S$empo_3)
s1<-filter(S,empo_3=='Water (non-saline)')
s2<-filter(S,empo_3=='Surface (saline)')
vs1<-dist(s1[,3:7],method='euclidean')
vs2<-dist(s2[,3:7],method='euclidean')
s1.dist<-dist(s1$sulfate_umol_per_l,method='euclidean')
s2.dist<-dist(s2$sulfate_umol_per_l,method='euclidean')
mantel(vs1,s1.dist,method='spearman')
mantel(vs2,s2.dist,method='spearman')



       
         
