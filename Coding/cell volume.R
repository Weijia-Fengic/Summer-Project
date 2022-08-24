#Data preparation 

##Data description
rm(list=ls())
setwd("F:/EA/Project/Madin")
sp_NCBI<-read.csv("F:/EA/Project/Madin/condensed_species_NCBI.csv")
str(sp_NCBI)
nrow(sp_NCBI)

#Select 'bacteria'
library("dplyr")
df<-filter(sp_NCBI,superkingdom=="Bacteria")
length(unique(df$phylum))
unique(df$cell_shape)
#create average cell size
df$d1_mid<-ifelse(!is.na(df$d1_up),(df$d1_lo+df$d1_up)/2,df$d1_lo)
df$d2_mid<-ifelse(!is.na(df$d2_up),(df$d2_lo+df$d2_up)/2,df$d2_lo)

#volume calculation
##Aggregate cell shapes for next steps (exclude 'NA'')
df$shapeagg<-df$cell_shape
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("coccus")] <- "spheroid"
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("bacillus","coccobacillus","vibrio")] <- "rod"
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("pleomorphic","filament", "star", "spiral", "irregular", "flask", "spindle", "fusiform", "disc", "disc ", "square", "branced", "triangular")] <- NA

#Create shapeagg column(include na)
df$shapeagg <- df$cell_shape
df$shapeagg[is.na(df$shapeagg)] <- "oval"
df$shapeagg[df$shapeagg %in% c("coccus")] <- "oval"
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("bacillus","coccobacillus","vibrio")] <- "rod"
df$shapeagg[!is.na(df$shapeagg) & df$shapeagg %in% c("pleomorphic","filament", "star", "spiral", "irregular", "flask", "spindle", "fusiform", "disc", "disc ", "square", "branced", "triangular")] <- NA


#calculate volumes
tmp <- df[!is.na(df$d1_mid) & !is.na(df$shapeagg) & df$shapeagg %in% c("rod","spheroid"), c("species_tax_id","d1_mid","d2_mid","shapeagg","cell_shape")]
tmp$volume <- NA

for (i in 1:nrow(tmp)) {
  if(!is.na(tmp$d1_mid[i])){
    if(tmp$shapeagg[i]=="spheroid") {
      
      d<-NA
    #If there is a difference between width (diameter) and length, use mean
      if(!is.na(tmp$d2_mid[i])) {
        #Calculate mean of the two values
        d <- (tmp$d1_mid[i]+tmp$d2_mid[i])/2
      } else {
        d <- tmp$d1_mid[i]
      }
      
      tmp$volume[i] <- 4/3*pi*(d/2)^3
      
    } else {
      if(!is.na(tmp$d2_mid[i])) {
        #Calculate as rod with hemispherical ends 
        
        #end volume:
        ends <- 4/3*pi*(tmp$d1_mid[i]/2)^3
        #body: length minus diameter used in ends
        body <- pi*(tmp$d1_mid[i]/2)^2*(tmp$d2_mid[i]-tmp$d1_mid[i])
        
        if(body>0) {
          tmp$volume[i] <- ends+body
        } else {
          tmp$volume[i] <- ends
        }
        
      }
    }
  }
}


tmp <- tmp[!is.na(tmp$volume),]
length(unique(tmp$species_tax_id))
log10(max(tmp$volume)*(10^(-18)))
log10(min(tmp$volume)*(10^(-18)))


##Get the value needed(include na)
tmp1 <- df[!is.na(df$d1_mid) & !is.na(df$shapeagg) & df$shapeagg %in% c("rod","oval"), c("species_tax_id","d1_mid","d2_mid","shapeagg","cell_shape")]
tmp1$volume <- NA
for(i in 1:nrow(tmp1)) {
  if(!is.na(tmp1$d1_mid[i])) {
    if(tmp1$shapeagg[i] == "oval") {
      d <- NA
      if(!is.na(tmp1$d2_mid[i])) {
        d1 <- tmp1$d1_mid[i]
        d2<-tmp1$d2_mid[i]
      } else {
        d1 <- tmp1$d1_mid[i]
        d2 <- tmp1$d1_mid[i]
      }
      tmp1$volume[i] <- 4/3*pi*(d1/2)^2*d2*(10^(-18))
    } else {
      if(!is.na(tmp1$d2_mid[i])) {
        #Calculate as rod with hemispherical ends
        #end volume:
        ends <- 4/3*pi*(tmp1$d1_mid[i]/2)^3*(10^(-18))
        #body: length minus diameter used in ends
        body <- pi*(tmp1$d1_mid[i]/2)^2*(tmp1$d2_mid[i]-tmp1$d1_mid[i])*(10^(-18))
        if(body>0) {
          tmp1$volume[i] <- ends+body
        } else {
          tmp1$volume[i] <- ends
        }
      }
    }
  }
}
tmp1 <- tmp1[!is.na(tmp1$volume),]
str(tmp1)

#Merge volume to main data frame
df <- df %>% left_join(tmp[,c("species_tax_id","volume")], by = "species_tax_id")
df1<-tmp1%>%left_join(df[,c("species_tax_id","isolation_source","species","genus","family","order","class","phylum")],by="species_tax_id")
nrow(df)
#volume distribution    
library("ggplot2")
vp<-ggplot(df,aes(x=log10(volume)))+
  geom_density(color="black",alpha=0.5)+
  xlab("Volume")+
  theme_classic()
ggsave(plot = vp, width = 8, height = 3, dpi = 300, filename = "volume distribution(Madin et al.).png")


#environment
env<- read.csv("F:/EA/Project/Madin/bacteria-archaea-traits-1.0.0/bacteria-archaea-traits-1.0.0/data/conversion_tables/environments.csv", as.is=TRUE)
df1<-df1%>%left_join(env[,c("Type","Cobo.Simon.habitat")],by = c("isolation_source"="Type"))
unique(df1$isolation_source)

#Create new habitat category 
##First letter upper case for environment
install.packages('stringr')
library('stringr')
df1$Cobo.Simon.habitat <- str_to_title(df1$Cobo.Simon.habitat, locale = "en")
#Use Cobo Simon habitat scheme
df1$habitat <- as.character(df1$Cobo.Simon.habitat)
unique(df1$habitat)

#Set any isolation_source that does not fit into the Cobo Simon scheme to "Other"
df1$habitat[!is.na(df1$isolation_source) & is.na(df1$habitat)] <- "Other"
df1$habitat[df1$habitat == "Therm"] <- "Thermal"
df1$habitat[grepl("plant|fungus|algae",df1$isolation_source)] <- "Other"
df1$habitat[grepl("feces|endotherm_surface|ectotherm_surface",df1$isolation_source)] <- "Other"
df1$habitat[df1$isolation_source %in% c("host","host_animal")] <- "Other"
#Where habitat is not "Other" split host into 'endo' and 'ecto'
df1$habitat[grepl("endotherm",df1$isolation_source) & !(df1$habitat == "Other")] <- "Endotherm"
df1$habitat[grepl("ectotherm",df1$isolation_source) & !(df1$habitat == "Other")] <- "Ectotherm"

#Set all species with no isolation_source/habitat information to habitat = "Other"
df1$habitat[is.na(df1$habitat)] <- "Other"
unique(df1$habitat)
df1$habitat[df1$habitat=="Low"]<-"Other"
unique(df1$habitat)
write.csv(df1,file = "volume_NCBI.csv")

#habitat and cell volume 
#Relationship in habitat and cell volume
ggplot(df,aes(x=habitat,y=log10(volume),fill=habitat))+
  geom_boxplot(color="black",alpha=0.5)+
  xlab("Habitat")+
  ylab("Log 10 volume")
  theme_classic()

#cell volume distribution (group by 'Habitat')   
ggplot(df,aes(x=habitat,y=log10(volume),fill=habitat))+
    geom_boxplot(color="black",alpha=0.5)+
    xlab("Habitat")+
    ylab("Log 10 volume")
  theme_classic()+
    theme(legend.text = element_text(size=25))+
    theme(legend.title = element_text(size=30))+
    theme(legend.position = c(0.7,0.5))+
    theme(axis.title = element_text(size=30))

