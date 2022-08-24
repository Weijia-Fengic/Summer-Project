#Cell size distribution 
rm(list=ls())
setwd("F:/EA/Project/EMP/Machado et al")
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")
library("ggpubr")
library("lme4")
library("gridExtra")
library("broom")
library("moments")


data<-read.csv("F:/EA/Project/EMP/Machado et al/EMP_Madin_gs.csv")
data1<-data%>%mutate(gs_Madin=coalesce(gs_Madin,gs_Machado))

#Genome size and cell volume of distribution in habitat
##volume distribution
habitat<-split(data,data$empo_3)
plot_list<-list()
pdf("v_empo3.pdf")
for (i in 1:length(habitat)) {
  sub_data1<-habitat[[i]]
  n1<-names(habitat[i])
  p1<-ggplot(sub_data1,aes(x=log10(volume),fill=empo_3))+
    geom_density(alpha=0.4)+
    xlab("Log 10 volume")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=22))+
    theme_classic()+
    labs(title=names(habitat[i]))
  plot_list[[i]]<-p1 
  print(plot_list[[i]])
  
}
dev.off()

##genome size distribution
plot_list1<-list()
pdf("gs_empo3.pdf")
for (i in 1:length(habitat)) {
  sub_data2<-habitat[[i]]
  n2<-names(habitat[i])
  p2<-ggplot(sub_data2,aes(x=gs_Madin/1000000,fill=empo_3))+
    geom_density(alpha=0.4)+
    xlab("Genome size (Mbp)")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=22))+
    theme_classic()+
    facet_wrap(.~empo_3,scales="free",ncol=2)+
    labs(title=names(habitat[i]))
  plot_list1[[i]]<-p2 
  print(plot_list1[[i]])
  
}
dev.off()

#Cell distribution within community (based on habitat)
##volume distribution (phylum by phylum)
phy1<-data1%>%filter(phylum=='Actinobacteria')
phy2<-data1%>%filter(phylum=='Bacteroidetes')
phy3<-data1%>%filter(phylum=='Firmicutes')
phy4<-data1%>%filter(phylum=='Proteobacteria')
p1<-ggplot(phy1,aes(x=log10(volume),fill=empo_3))+
    geom_density(alpha=0.8)+
    xlab("Log 10 volume")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=22))+
    theme_classic()+
    labs(title='Actinobacteria',
         fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

p2<-ggplot(phy2,aes(x=log10(volume),fill=empo_3))+
  geom_density(alpha=0.8)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()+
  labs(title='Bacteroidetes',
       fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

p3<-ggplot(phy3,aes(x=log10(volume),fill=empo_3))+
  geom_density(alpha=0.8)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()+
  labs(title='Firmicutes',
       fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

p4<-ggplot(phy4,aes(x=log10(volume),fill=empo_3))+
  geom_density(alpha=0.8)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()+
  labs(title='Proteobacteria',
       fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

ggarrange(p1,p2,p3,p4,ncol = 2,nrow=2,labels = c("a","b","c","d"), common.legend = T,legend = "right")

#genome size distribuion
p5<-ggplot(phy1,aes(x=log10(gs_Madin),fill=empo_3))+
  geom_density(alpha=0.8)+
  xlab("Log 10 genome size")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()+
  labs(title='Actinobacteria',
       fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

p6<-ggplot(phy2,aes(x=log10(gs_Madin),fill=empo_3))+
  geom_density(alpha=0.8)+
  xlab("Log 10 genome size")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()+
  labs(title='Bacteroidetes',
       fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

p7<-ggplot(phy3,aes(x=log10(gs_Madin),fill=empo_3))+
  geom_density(alpha=0.8)+
  xlab("Log 10 genome size")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()+
  labs(title='Firmicutes',
       fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

p8<-ggplot(phy4,aes(x=log10(gs_Madin),fill=empo_3))+
  geom_density(alpha=0.8)+
  xlab("Log 10 genome size")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()+
  labs(title='Proteobacteria',
       fill='Habitat')+
  scale_fill_manual(values=getPalette1(colourCount1))

ggarrange(p5,p6,p7,p8,ncol = 2,nrow=2,labels = c("a","b","c","d"), common.legend = T,legend = "right")


#Model fitting: habitat level community
result<-lapply(habitat,function(data) lm(volume~gs_Madin,data = data))
lapply(result,summary)
#Extract the values
data$x<-data$gs_Madin
data$y<-data$volume
Coef1<-data %>% group_by(empo_3) %>% 
  summarise(
    r.sqr = summary(lm(y~x))$r.squared,
    Intercept=summary(lm(y~x))$coefficients[[1]],
    slop = summary(lm(y~x))$coefficients[[2]],
    p_value=summary(lm(y~x))$coefficients[[3]])
newdata1<-merge(Coef1,data,by="empo_3")
newdata2<-newdata1%>%
  group_by(empo_3)%>%
  summarise(R2=sum(r.sqr),
            intercept=sum(Intercept),
            Slop=sum(slop),
            Pr=sum(p_value))

habitatcom<-write.table (newdata2, file ="habitat community model summary.csv", sep =",")

#genome size and volume relationship in each habitat
plot_list2<-list()
pdf("gs&V_empo3.pdf")
for (i in 1:length(habitat)) {
  sub_data<-habitat[[i]]
  n<-names(habitat[i])
  p<-ggplot(sub_data,aes(x=log10(volume),y=gs_Madin/1000000,color=phylum))+
    geom_point(size=3,alpha=0.4)+
    xlab("Log 10 volume")+
    ylab("Genome size (Mbp)")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=22))+
    theme_classic()+
    labs(title=names(habitat[i]))
  plot_list2[[i]]<-p 
  print(plot_list2[[i]])
  
}
dev.off()

#Density plots: single phylum group in each habitat
#volume distribution
length(unique(data$phylum))
phylum<-split(data,data$phylum)
plot_list<-list()
pdf("volume_phylum.pdf")
for (i in 1:length(phylum)) {
  sub_data<-phylum[[i]]
  n<-names(phylum[i])
  p1<-ggplot(sub_data,aes(x=log10(volume),fill=empo_3))+
    geom_density(alpha=0.4)+
    xlab("Log 10 volume")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=22))+
    theme_classic()+
    facet_wrap(.~empo_3,scales="free",ncol=3)+
    labs(title=names(phylum[i]))
  plot_list[[i]]<-p1 
  print(plot_list[[i]])
  
}
dev.off()

#genome size distribuion
plot_list<-list()
pdf("gs_phylum.pdf")
for (i in 1:length(phylum)) {
  sub_data<-phylum[[i]]
  n<-names(phylum[i])
  p1<-ggplot(sub_data,aes(x=gs_Madin/1000000,fill=empo_3))+
    geom_density(alpha=0.4)+
    xlab("Genome size (Mbp)")+
    theme(axis.title = element_text(size=25),
          axis.text = element_text(size=20),
          legend.text = element_text(size=18),
          legend.title = element_text(size=22))+
    theme_classic()+
    facet_wrap(.~empo_3,scales="free",ncol=3)+
    labs(title=names(phylum[i]))
  plot_list[[i]]<-p1 
  print(plot_list[[i]])
  
}
dev.off()

#Model fitting: sample level community
list_of_Sample<-split(data1,data1$sample)
results<-lapply(list_of_Sample,function(data1) lm(log10(gs_Madin)~log10(volume),data = data1))
lapply(results,summary)
#Extract the values and basic statistical description for sample community
data1$x<-log10(data1$volume)
data1$y<-log10(data1$gs_Madin)
Coef2<-data1 %>% group_by(sample) %>% 
  summarise(
    r.sqr = summary(lm(y~x))$r.squared,
    Intercept=summary(lm(y~x))$coefficients[1,1],
    slop = summary(lm(y~x))$coefficients[1,2],
    p_value=summary(lm(y~x))$coefficients[,"Pr(>|t|)"])
newdata3<-merge(Coef2,data1,by="sample")
newdata4<-newdata3%>%
  group_by(sample)%>%
  summarise(R2=r.sqr,
            intercept=Intercept,
            Slop=slop,
            Pr=p_value,
            GS_mean=mean(gs_Madin),
            GS_median=median(gs_Madin),
            GS_var=var(gs_Madin),
            GS_max=max(gs_Madin),
            GS_min=min(gs_Madin),
            GS_skew=skewness(gs_Madin),
            GS_kurt=kurtosis(gs_Madin),
            V_mean=mean(volume),
            V_median=median(volume),
            V_var=var(volume),
            V_max=max(volume),
            V_min=min(volume),
            V_skew=skewness(volume),
            V_kurt=kurtosis(volume),
            comm_size=length(unique(Species)))%>%
  distinct(sample,.keep_all = T) #remove duplicate rows across 'sample' column  
            
#merge metadata and sample community
metadata<-read_tsv("F:/EA/Project/EMP/Machado et al/data/data/emp_qiime_mapping_qc_filtered.tsv")
metadata<-metadata%>%rename("sample"="#SampleID")
newdata4<-newdata4%>%left_join(metadata[,c("sample","empo_3","ph","latitude_deg",
                                           "longitude_deg","depth_m","altitude_m", 
                                           "elevation_m","adiv_shannon","temperature_deg_c",
                                           "oxygen_mg_per_l","phosphate_umol_per_l","ammonium_umol_per_l",
                                           "nitrate_umol_per_l","sulfate_umol_per_l")],by="sample")

samplecom<-write.table (newdata4, file ="sample_community_summary.csv", sep =",")

#volume and genome size distribution in sample community (examples)
sample1<-data%>%filter(sample=='2192.H03a.Hand.1519.lane6.NoIndex.L006')
s1_GS1<-ggplot(sample1,aes(x=gs_Madin/1000000))+
  geom_density(alpha=0.4,size=2)+
  xlab("Genome size (Mbp)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()
s1_V1<-ggplot(sample1,aes(x=log10(volume)))+
  geom_density(alpha=0.4,size=2)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()
  
sample2<-data%>%filter(sample=='905.PLO6.t2.OSCANOX.0to1cm.1')
s2_GS1<-ggplot(sample2,aes(x=gs_Madin/1000000))+
  geom_density(alpha=0.4,size=2)+
  xlab("Genome size (Mbp)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()
s2_V1<-ggplot(sample2,aes(x=log10(volume)))+
  geom_density(alpha=0.4,size=2)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()

sample3<-data%>%filter(sample=='945.P11.B4.lane3.NoIndex.L003')
s3_GS1<-ggplot(sample3,aes(x=gs_Madin/1000000))+
  geom_density(alpha=0.4,size=2)+
  xlab("Genome size (Mbp)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()
s3_V1<-ggplot(sample3,aes(x=log10(volume)))+
  geom_density(alpha=0.4,size=2)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()

sample4<-data%>%filter(sample=='1747.DZF.6252012.A.shipping.box')
s4_GS1<-ggplot(sample4,aes(x=gs_Madin/1000000))+
  geom_density(alpha=0.4,size=2)+
  xlab("Genome size (Mbp)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()
s4_V1<-ggplot(sample4,aes(x=log10(volume)))+
  geom_density(alpha=0.4,size=2)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()

sample5<-data%>%filter(sample=='2382.GM.181.R5.bulk.10.12.lane7.NoIndex.L007.sequences')
s5_GS1<-ggplot(sample5,aes(x=gs_Madin/1000000))+
  geom_density(alpha=0.4,size=2)+
  xlab("Genome size (Mbp)")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()
s5_V1<-ggplot(sample5,aes(x=log10(volume)))+
  geom_density(alpha=0.4,size=2)+
  xlab("Log 10 volume")+
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=22))+
  theme_classic()




