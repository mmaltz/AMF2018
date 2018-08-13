#Create SSU boxplots of functional groups
#Author : Michala Phillips
library(ggplot2)
setwd("D:/R_Workspace/Amplicon_stuff/Chapter_2_NAUData/data/SSU_Redo/")
ssu<-read.csv("2018_08_02_ssu_funcguild_read_rich.csv")
names(ssu)[names(ssu) == "TREATMENT"] <- "Host.Plant.Type"

roots <- subset(ssu, ssu$TYPE=="ROOT")
soil <- subset(ssu, ssu$TYPE=="SOIL")


richbox <- ggplot(soil, aes(x = Functional.Group , y = OTU_Richness_Sample)) + geom_boxplot(aes(fill = Host.Plant.Type)) 
richbox <- richbox + theme_bw(base_size = 15) + xlab("Functional Group") +ylab("Soil AMF Taxa Richness")    
richbox <-richbox + scale_y_continuous(limits=c(0,150))
richbox <- richbox + scale_fill_manual(values=c("#999999", "#FFFFFF"))
richbox



#root richness
richbox <- ggplot(roots, aes(x = Functional.Group , y = OTU_Richness_Sample)) + geom_boxplot(aes(fill = Host.Plant.Type)) 
richbox <- richbox + theme_bw(base_size = 15)+ xlab("Functional Group") +ylab("Root AMF Taxa Richness")    
richbox <-richbox + scale_y_continuous(limits=c(0,150)) +annotate("text", x = 3, y = 75, label = "**", size=10)
richbox <- richbox + scale_fill_manual(values=c("#999999", "#FFFFFF"))
richbox 


#soil read
readbox <- ggplot(soil, aes(x = Functional.Group , y = Read_Abundance_Sample)) + geom_boxplot(aes(fill = Host.Plant.Type)) 
readbox <- readbox + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line.x = element_line(color="black", size = 1),
                           axis.line.y = element_line(color="black", size = 1))+ xlab("Functional Group") +ylab("AMF Read Abundance")    
readbox <-readbox + scale_y_continuous(limits=c(0,500))
readbox 

#root read
readbox <- ggplot(roots, aes(x = Functional.Group , y = Read_Abundance_Sample)) + geom_boxplot(aes(fill = Host.Plant.Type)) 
readbox <- readbox + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.background = element_blank(), axis.line.x = element_line(color="black", size = 1),
                           axis.line.y = element_line(color="black", size = 1))+ xlab("Functional Group") +ylab("AMF Read Abundance")    
readbox <-readbox + scale_y_continuous(limits=c(0,500))
readbox 




 