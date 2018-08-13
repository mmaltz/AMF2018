

#load Required packages
library(dplyr)
library(reshape2)
library(ggplot2)

#Set wd and import data files
setwd("D:/R_Workspace/Amplicon_stuff/Chapter_2_NAUData/data/SSU_redo")
ssu<-read.csv("SSU_Redo_CSS_table_sorted_with_metadata.csv", header = TRUE)
metassu<-read.csv("MP_SSU_MAP.csv")

#Rename SSU hedings
colnames(ssu) <- c("ssuotu", 1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 20, 
                   25, 26, 27, 28, 33, 34, 35, 36, 41, 42, 43, 44, 49, 
                   50, 51, 52, 57, 58, 59, 60, 65, 66, 67, 68, 74, 75, 
                   76, 81, 82, 83, 84, 89, 90, 92,"ssutaxonomy","ssukingdom", 
                   "ssuphylum", "ssuclass", "ssuorder", "ssufamily", "ssugenus", "ssuspecies" )
#remove Geosiphonaceae
drop.family <- c(
  "Geosiphonaceae"
)
ssu  <- ssu[-which(ssu$ssufamily %in% drop.family),]

drop <-c("No blast hit")
ssu  <- ssu[-which(ssu$ssutaxonomy %in% drop),]


#add in fucntional groups
for (i in 1:179){
  
  if (ssu$ssufamily[i]=="Glomeraceae") {
    ssu$Functional.Group[i]<-"Rhizophilic"
  } else if ( ssu$ssufamily[i]=="Claroideoglomeraceae") {
    ssu$Functional.Group[i]<-"Rhizophilic"
  } else if(ssu$ssufamily[i]=="Gigasporaceae"){
    ssu$Functional.Group[i]<-"Edaphophilic"
  } else if(ssu$ssufamily[i]=="Diversisporaceae") {
    ssu$Functional.Group[i]<-"Edaphophilic"
  } else if (ssu$ssufamily[i]=="Paraglomeraceae") {
    ssu$Functional.Group[i]<-"Rhizophilic"
  } else if(ssu$ssufamily[i]=="Ambisporaceae"){
    ssu$Functional.Group[i]<-"Ancestral"
  } else if(ssu$ssufamily[i]=="Archaeosporaceae"){
    ssu$Functional.Group[i]<-"Ancestral"
  } else if(ssu$ssufamily[i]=="Acaulosporaceae"){
    ssu$Functional.Group[i]<-"Ancestral"
  } 
  
}

#Take data from wide to long and renaming variables

ssul <-melt(ssu)
names(ssul)[names(ssul)=="variable"] <- "ID1"
names(ssul)[names(ssul)=="value"] <- "ssureads"


#Functional group OTU richnesses and read abundance for each functional 
#group I LIKE THIS WAY A LOT BECAUSE YOU CAN NAME THE ACTUAL HEADERS
#Below is in format for graphing
ssuguildREADRICH<- data.frame(ssul %>%
                                group_by(ID1,Functional.Group) %>%
                                summarise(OTU_Richness_Sample = length(unique(ssuotu[ssureads>0])),
                                          Read_Abundance_Sample = sum(ssureads)
                                ) 
                              
)




#Make SSU data frame with family level OTU richness and Read abundance for each family

ssufamily <- data.frame(ssul %>%
                          group_by(ID1, ssufamily) %>%
                          summarise(Family_OTU_Richness = length(unique(ssuotu[ssureads>0])), 
                                    Family_Read_Abundance = sum(ssureads)
                          ))

#Add metadata into ssu family

names(metassu)[names(metassu)=="ï..SampleID"] <- "ID1"
ssufamily<-merge(ssufamily,metassu, by="ID1")
ssul<-merge(ssul,metassu,by="ID1")

#Adds meta data to data for box plot


ssuguildREADRICH<-merge(ssuguildREADRICH,metassu, by="ID1")
#write.csv(ssuguildREADRICH, file = "2018_08_02_ssu_funcguild_read_rich.csv")
###For ssufamily, moves ID to single line and family to the length
LssufamilyOTURICH <- dcast(ssufamily, ID ~ ssufamily, value.var = "Family_OTU_Richness",  fun.aggregate = sum )
LssufamilyOTUREAD <- dcast(ssufamily, ID ~ ssufamily, value.var = "Family_Read_Abundance",  fun.aggregate = sum )

#write.csv(LssufamilyOTURICH, file ="2018_08_02_ssufamily_rich.csv")
#write.csv(LssufamilyOTUREAD, file ="2018_08_02_ssufamily_read.csv")

#Functional group read and richness together

ssuguildREADRICHlong<- data.frame(ssul %>%
                                    group_by(ID) %>%
                                    summarise(OTU_Richness_Sample = length(unique(ssuotu[ssureads>0])),
                                              Ancestral_Richness = length(unique(ssuotu[ssureads>0 & Functional.Group =="Ancestral"])),
                                              Edaphophilic_Richness = length(unique(ssuotu[ssureads>0 & Functional.Group =="Edaphophilic"])),
                                              Rhizophilic_Richness = length(unique(ssuotu[ssureads>0 & Functional.Group =="Rhizophilic"])),
                                              Ancestral_Read_Abundance = sum(ssureads[Functional.Group == "Ancestral"]),
                                              Edaphophilic_Read_Abundance = sum(ssureads[Functional.Group == "Edaphophilic"]),
                                              Rhizophilic_Read_Abundance = sum(ssureads[Functional.Group == "Rhizophilic"]),
                                              Read_Abundance_Sample = sum(ssureads)
                                    ) 
)

ssuguildREADRICHlong <- merge(ssuguildREADRICHlong, metassu, by = "ID")                             
#write.csv(ssuguildREADRICHlong, file = "2018_08_02_ssu_funcguild_read_rich.csv")

#Functional group read and richness separate 

ssuguildREAD<- data.frame(ssul %>%
                            group_by(ID) %>%
                            summarise(Ancestral_Read_Abundance = sum(ssureads[Functional.Group == "Ancestral"]),
                                      Edaphophilic_Read_Abundance = sum(ssureads[Functional.Group == "Edaphophilic"]),
                                      Rhizophilic_Read_Abundance = sum(ssureads[Functional.Group == "Rhizophilic"])
                            ) 
)

ssuguildRICH<- data.frame(ssul %>%
                            group_by(ID) %>%
                            summarise(OTU_Richness_Sample = length(unique(ssuotu[ssureads>0])),
                                      Ancestral_Richness = length(unique(ssuotu[ssureads>0 & Functional.Group =="Ancestral"])),
                                      Edaphophilic_Richness = length(unique(ssuotu[ssureads>0 & Functional.Group =="Edaphophilic"])),
                                      Rhizophilic_Richness = length(unique(ssuotu[ssureads>0 & Functional.Group =="Rhizophilic"]))
                            ) 
)

#Add metadata to above 

ssuguildRICH<-merge(ssuguildRICH, metassu, by="ID")
ssuguildREAD<-merge(ssuguildREAD, metassu, by="ID")


#Separating data by root and soil

ssuroot.readreach<-subset(ssuguildREADRICH,ssuguildREADRICH$TYPE == "ROOT")
ssuguildREADRICHlong<-merge(ssuguildREADRICHlong,metassu,by="ID")
ssuroot.readreachlong<-subset(ssuguildREADRICHlong,ssuguildREADRICHlong$TYPE == "ROOT")
#soil
ssusoil.readreach<-subset(ssuguildREADRICH,ssuguildREADRICH$TYPE == "SOIL")
ssuguildREADRICHlong<-merge(ssuguildREADRICHlong,metassu,by="ID")
ssusoil.readreachlong<-subset(ssuguildREADRICHlong,ssuguildREADRICHlong$TYPE == "SOIL")

#Pull out roots and soil for Family DFs

ssufamroot<-subset.data.frame(ssufamily,TYPE == "ROOT")
names(ssufamroot)[names(ssufamroot)=="ssufamily"]<-"Family"
ssufamsoil<-subset.data.frame(ssufamily,TYPE == "SOIL")
names(ssufamsoil)[names(ssufamsoil)=="ssufamily"]<-"Family"




