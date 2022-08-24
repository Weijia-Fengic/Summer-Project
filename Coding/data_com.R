#Data combination (Machado et al., Zhaoyi, Madin et. al)
rm(list=ls())
setwd("F:/EA/Project/EMP/Machado et al")
library("readr")
library("dplyr")
library("tidyr")
library("tidyverse")
library("ggplot2")

#Load sample data and metadata
s<-read_tsv("F:/EA/Project/EMP/Machado et al/data/data/emp_150bp_filtered.tsv")
s['value']<-1
s<-s%>%separate(col=org_id,into = c("genus","species","ID"),sep="_")
metadata<-read_tsv("F:/EA/Project/EMP/Machado et al/data/data/emp_qiime_mapping_qc_filtered.tsv")
metadata<-metadata%>%rename("sample"="#SampleID")

#merge sample data and metadata
s<-s%>%left_join(metadata[,c("sample","empo_1","empo_2","empo_3","ph")],by="sample")

#merge genome size(Madin, Machado) data and sample data
#merge genome size from Madin et al at first
sp_NCBI<-read.csv("F:/EA/Project/Madin/condensed_species_NCBI.csv")
df<-filter(sp_NCBI,superkingdom=="Bacteria")
df1<-df[,c("species","genome_size")]
df1<-df1%>%separate(col=species,into = c("genus","species"),sep = " ")
genome<-merge(s,df1,by.x = c("genus","species"),by.y = c("genus","species"))
genome<-genome%>%rename("gs_Madin"="genome_size")

#merge genome size from Machado et al.
gM<-read_tsv("F:/EA/Project/EMP/Machado et al/data/data/genome_sizes.tsv")
gM1<-gM[,c("bp","org_id")]
gM1<-gM1%>%separate(col=org_id,into = c("genus","species","ID"),sep="_")
genome<-genome%>%left_join(gM1[,c("genus","species","bp")],by=c("genus","species"))
genome<-genome%>%rename("gs_Machado"="bp")

#merge volume and sample data
#Volume from Madin et al and sample data
MV<-read.csv("F:/EA/Project/Madin/volume_NCBI.csv")  #Only madin et al. data
subMV<-MV[,c("species","shapeagg","cell_shape","volume","phylum")]
subMV<-subMV%>%separate(col=species,into = c("genus","species"),sep=" ")
data<-merge(genome,subMV,by.x = c("genus","species"),by.y = c("genus","species"))
data$Species<-paste(data$genus,data$species,sep=" ")
data<-data%>%left_join(metadata[,c("sample","latitude_deg",
                                  "longitude_deg","depth_m","altitude_m", 
                                  "elevation_m","adiv_shannon","temperature_deg_c",
                                  "oxygen_mg_per_l","phosphate_umol_per_l","ammonium_umol_per_l",
                                  "nitrate_umol_per_l","sulfate_umol_per_l")],by="sample")
write.csv(data,"EMP_Madin_gs.csv")
#Basic data description 
#No of sample
length(unique(data$sample))
#No of species
length(unique(data$Species))

#Genome size density plot
d<-ggplot(data)+
  geom_density(aes(x=gs_Madin),color="black",fill="#FF0000",alpha=0.4,size=1,show.legend = T)+
  geom_density(aes(x=gs_Machado),color="black",fill="#36BED9",alpha=0.4,size=1,show.legend = T)+
  xlab("Genome size (Mbp)")+
  theme_classic()+
  scale_color_manual(name=NULL,
                     guide="legend",
                    breaks =c("Madin","Machado"),
                     values=c("Madin"="#FF0000","Machado"="#36BED9"))+
 theme(axis.title = element_text(size=27),
       axis.text = element_text(size=20),
       legend.text = element_text(size=25),
       legend.position = c(0.7,0.6))


#Volume and environment 
ggplot(data,aes(empo_3,log10(volume)))+
  geom_boxplot(aes(fill=empo_3))+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90),
        axis.title = element_text(size=23),
        axis.text = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=12))

#Compare genome size from Machado and Madin
d1<-ggplot(data,aes(gs_Machado,gs_Madin))+
  geom_point(fill="black")+
  theme_classic()+
  theme(axis.text.x=element_text(angle = 90),
        axis.title = element_text(size=23),
        axis.text = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=12))

                               
                              
 

