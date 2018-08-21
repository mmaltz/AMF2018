 #GLMM for roots and soil OTU Richness SSU by functional group-modified from "GLM_NAU_6_2017"
#---different models for each functional group
#Author: Michala Phillips

#load packages
library(car)
require(MASS)
library(lme4)

setwd("D:/R_Workspace/Amplicon_stuff/Chapter_2_NAUData/data/SSU_redo")
ssu<-read.csv("2018_08_02_ssu_funcguild_read_rich.csv")
phnut<-read.csv("Nutrient_ph_nausamps1.csv", header = T)
#add nutrient data etc
ssu<-merge(ssu, phnut, by = "ID")

#ssu$Rhizophilic <- ifelse(ssu$Functional.Group == "Rhizophilic", 1, 0)
#ssu$Late <- ifelse(ssu$Functional.Group == "Late", 1, 0)
#ssu$Mixed <- ifelse(ssu$Functional.Group == "Mixed", 1, 0)


roots<-subset(ssu, ssu$TYPE=="ROOT")

#subsets roots by functional group
rootsear <- subset(roots, roots$Functional.Group =="Rhizophilic")
rootsmixed <-  subset(roots, roots$Functional.Group == "Ancestral")
rootslate <- subset(roots, roots$Functional.Group == "Edaphophilic")
##soilsubset
soil<-subset(ssu, ssu$TYPE=="SOIL")

#subsets soil by functional group
soilear <- subset(soil, soil$Functional.Group =="Rhizophilic")
soilmixed <-  subset(soil, soil$Functional.Group == "Ancestral")
soillate <- subset(soil, soil$Functional.Group == "Edaphophilic")
#############OTU RICHNESS
##############EARLY####################

# find the right distribution for early roots
rootsear$OTU_Richness_Sample.t <- rootsear$OTU_Richness_Sample + 1
qqp(rootsear$OTU_Richness_Sample.t, "norm")

#LOGNORMAL
qqp(rootsear$OTU_Richness_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(rootsear$OTU_Richness_Sample.t, "Negative Binomial")
qqp(rootsear$OTU_Richness_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(rootsear$OTU_Richness_Sample.t, "Poisson")
qqp(rootsear$OTU_Richness_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(rootsear$OTU_Richness_Sample.t, "gamma")
qqp(rootsear$OTU_Richness_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])
#
earlyrich <- glm(OTU_Richness_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                 + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3 + TREATMENT*pH 
                 + TREATMENT*P + SITE*pH, data = rootsear, family = gaussian(link = "identity"))
stepAIC(earlyrich)
earlyrichfinal <- glm(OTU_Richness_Sample.t ~ TREATMENT + SITE + P + 
                        NH4 + NO3 + pH + TREATMENT:SITE + SITE:NH4 + TREATMENT:NH4,
                      family = gaussian(link = "identity")
                      , data = rootsear)

summary(earlyrichfinal)
plot(earlyrichfinal)


#soil early

# find the right distribution for early soil
soilear$OTU_Richness_Sample.t <- soilear$OTU_Richness_Sample + 1
qqp(soilear$OTU_Richness_Sample.t, "norm")

#LOGNORMAL
qqp(soilear$OTU_Richness_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(soilear$OTU_Richness_Sample.t, "Negative Binomial")
qqp(soilear$OTU_Richness_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(soilear$OTU_Richness_Sample.t, "Poisson")
qqp(soilear$OTU_Richness_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(soilear$OTU_Richness_Sample.t, "gamma")
qqp(soilear$OTU_Richness_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#check that mean is above 5 for PQL GLM
mean(soilear$OTU_Richness_Sample)

th<-theta.ml(y =soilear$OTU_Richness_Sample.t, mu = mean(soilear$OTU_Richness_Sample.t, limit = 10))

earlyrich <- glm(OTU_Richness_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                 + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3
                 ,data = soilear,family = gaussian(link = "identity"))
stepAIC(earlyrich)
earlyrichfinal <- glm(OTU_Richness_Sample.t ~ TREATMENT + SITE + P + NH4 + SITE:P, data = soilear, family = gaussian(link = "identity"))
summary(earlyrichfinal)
plot(earlyrichfinal)
########################## Mixed ###########################
# find the right distribution for mixed roots
rootsmixed$OTU_Richness_Sample.t <- rootsmixed$OTU_Richness_Sample + 1
qqp(rootsmixed$OTU_Richness_Sample.t, "norm")

#LOGNORMAL
qqp(rootsmixed$OTU_Richness_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(rootsmixed$OTU_Richness_Sample.t, "Negative Binomial")
qqp(rootsmixed$OTU_Richness_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(rootsmixed$OTU_Richness_Sample.t, "Poisson")
qqp(rootsmixed$OTU_Richness_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(rootsmixed$OTU_Richness_Sample.t, "gamma")
qqp(rootsmixed$OTU_Richness_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


mixedrich <- glm(OTU_Richness_Sample.t ~TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                 + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3
                 ,data = rootsmixed,family = gaussian(link = "identity"))
stepAIC(mixedrich)
mixedrichfinal <- glm(OTU_Richness_Sample.t ~ TREATMENT + SITE + P + NH4 + NO3 + pH + TREATMENT:NH4 + TREATMENT:NO3 
                      ,data = rootsmixed,family = gaussian(link = "identity"))
summary(mixedrichfinal)
plot(mixedrichfinal)



# find the right distribution for mixedly soil
soilmixed$OTU_Richness_Sample.t <- soilmixed$OTU_Richness_Sample + 1
qqp(soilmixed$OTU_Richness_Sample.t, "norm")

#LOGNORMAL
qqp(soilmixed$OTU_Richness_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(soilmixed$OTU_Richness_Sample.t, "Negative Binomial")
qqp(soilmixed$OTU_Richness_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(soilmixed$OTU_Richness_Sample.t, "Poisson")
qqp(soilmixed$OTU_Richness_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(soilmixed$OTU_Richness_Sample.t, "gamma")
qqp(soilmixed$OTU_Richness_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#check that mean is above 5 for PQL GLM
mean(soilmixed$OTU_Richness_Sample)

mixedlyrich <- glm(OTU_Richness_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                   + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3, data = soilmixed, family = gaussian(link ="identity"))

stepAIC(mixedlyrich)
mixedlyrichfinal <- glm(OTU_Richness_Sample.t ~ SITE + P + NH4 + NO3 + SITE:NH4 + SITE:P, family = gaussian(link = "identity"), data = soilmixed)


summary(mixedlyrichfinal)
plot(mixedlyrichfinal)

########################## Late ###################

# find the right distribution for late roots
rootslate$OTU_Richness_Sample.t <- rootslate$OTU_Richness_Sample + 1
qqp(rootslate$OTU_Richness_Sample.t, "norm")

#LOGNORMAL
qqp(rootslate$OTU_Richness_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(rootslate$OTU_Richness_Sample.t, "Negative Binomial")
qqp(rootslate$OTU_Richness_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(rootslate$OTU_Richness_Sample.t, "Poisson")
qqp(rootslate$OTU_Richness_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(rootslate$OTU_Richness_Sample.t, "gamma")
qqp(rootslate$OTU_Richness_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

th<-theta.ml(y =rootslate$OTU_Richness_Sample.t, mu = mean(rootslate$OTU_Richness_Sample.t, limit = 10))

laterich <- glm(OTU_Richness_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3 + TREATMENT*pH 
                + TREATMENT*P + SITE*pH, 
                family =  neg.bin(theta = th), data = rootslate)

stepAIC(laterich)
laterichfinal <- glm(OTU_Richness_Sample.t ~ TREATMENT + SITE + NH4 + 
                       NO3 + pH + SITE:NH4 + TREATMENT:NH4 + SITE:pH, family = neg.bin(theta = th), 
                     data = rootslate)
summary(laterichfinal)
plot(laterichfinal)



# find the right distribution for late soil
soillate$OTU_Richness_Sample.t <- soillate$OTU_Richness_Sample + 1
qqp(soillate$OTU_Richness_Sample.t, "norm")

#LOGNORMAL
qqp(soillate$OTU_Richness_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(soillate$OTU_Richness_Sample.t, "Negative Binomial")
qqp(soillate$OTU_Richness_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(soillate$OTU_Richness_Sample.t, "Poisson")
qqp(soillate$OTU_Richness_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(soillate$OTU_Richness_Sample.t, "gamma")
qqp(soillate$OTU_Richness_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

latelyrich <- glm(OTU_Richness_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                  + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3, data = soillate, family = gaussian(link ="identity"))
stepAIC(latelyrich)
laterichfinal <- glm(OTU_Richness_Sample.t ~  TREATMENT + SITE + P + NH4 + NO3 + TREATMENT:SITE + 
                       SITE:NH4 + SITE:P + SITE:NO3 + TREATMENT:NO3
                     , data =soillate, family = gaussian(link ="identity"))
plot(laterichfinal)
summary(laterichfinal)

################READ ABUNDANCE#############
##############EARLY####################
# find the right distribution for early roots
rootsear$Read_Abundance_Sample.t <- rootsear$Read_Abundance_Sample + 1
qqp(rootsear$Read_Abundance_Sample.t, "norm")

#LOGNORMAL
qqp(rootsear$Read_Abundance_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(rootsear$Read_Abundance_Sample.t, "Negative Binomial")
qqp(rootsear$Read_Abundance_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(rootsear$Read_Abundance_Sample.t, "Poisson")
qqp(rootsear$Read_Abundance_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(rootsear$Read_Abundance_Sample.t, "gamma")
qqp(rootsear$Read_Abundance_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#
earlyrich <- glm(Read_Abundance_Sample.t ~TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                 + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3 + TREATMENT*pH 
                 + TREATMENT*P + SITE*pH, data = rootsear, 
                 family = gaussian(link = "identity"))
stepAIC(earlyrich)
earlyrichfinal <- glm(Read_Abundance_Sample.t ~TREATMENT + SITE + P + NH4 + NO3 + 
                        pH + TREATMENT:SITE + SITE:NH4 + TREATMENT:NH4, data = rootsear, 
                      family = gaussian(link = "identity"))
#family = neg.bin(theta = th))
summary(earlyrichfinal)
plot(earlyrichfinal)

#soil early

# find the right distribution for early soil
soilear$Read_Abundance_Sample.t <- soilear$Read_Abundance_Sample + 1
qqp(soilear$Read_Abundance_Sample.t, "norm")

#LOGNORMAL
qqp(soilear$Read_Abundance_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(soilear$Read_Abundance_Sample.t, "Negative Binomial")
qqp(soilear$Read_Abundance_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(soilear$Read_Abundance_Sample.t, "Poisson")
qqp(soilear$Read_Abundance_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(soilear$Read_Abundance_Sample.t, "gamma")
qqp(soilear$Read_Abundance_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


earlyrich <- glm(Read_Abundance_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                 + SITE*P + TREATMENT*NH4,family = gaussian(link = "identity"), data = soilear)
stepAIC(earlyrich)
earlyrich <- glm(Read_Abundance_Sample.t ~TREATMENT + SITE + NH4 + 
                   pH + TREATMENT:SITE + SITE:NH4 ,family = gaussian(link = "identity"), data = soilear)
summary(earlyrich)
plot(earlyrich)
########################## Mixed ###########################
# find the right distribution for mixed roots
rootsmixed$Read_Abundance_Sample.t <- rootsmixed$Read_Abundance_Sample + 1
qqp(rootsmixed$Read_Abundance_Sample.t, "norm")

#LOGNORMAL
qqp(rootsmixed$Read_Abundance_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(rootsmixed$Read_Abundance_Sample.t, "Negative Binomial")
qqp(rootsmixed$Read_Abundance_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(rootsmixed$Read_Abundance_Sample.t, "Poisson")
qqp(rootsmixed$Read_Abundance_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(rootsmixed$Read_Abundance_Sample.t, "gamma")
qqp(rootsmixed$Read_Abundance_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#check that mean is above 5 for PQL GLM
mean(rootsmixed$Read_Abundance_Sample)

th<-theta.ml(y =rootsmixed$Read_Abundance_Sample.t, mu = mean(rootsmixed$Read_Abundance_Sample.t, limit = 10))

mixedrich <- glm(Read_Abundance_Sample.t ~TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                 + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3 + TREATMENT*pH 
                 + TREATMENT*P + SITE*pH, data = rootsmixed, family = gaussian(link = "identity"))
#, 
#family = gaussian(link = "log"))
stepAIC(mixedrich)
mixedrichfinal <- glm(Read_Abundance_Sample.t ~  SITE +  TREATMENT + P + 
                        NH4 + NO3 + pH + TREATMENT:SITE + SITE:NH4 + TREATMENT:P + 
                        SITE:pH, family = gaussian(link = "identity"), data = rootsmixed)
summary(mixedrichfinal)
plot(mixedrichfinal)

# find the right distribution for mixedly soil
soilmixed$Read_Abundance_Sample.t <- soilmixed$Read_Abundance_Sample + 1
qqp(soilmixed$Read_Abundance_Sample.t, "norm")

#LOGNORMAL
qqp(soilmixed$Read_Abundance_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(soilmixed$Read_Abundance_Sample.t, "Negative Binomial")
qqp(soilmixed$Read_Abundance_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(soilmixed$Read_Abundance_Sample.t, "Poisson")
qqp(soilmixed$Read_Abundance_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(soilmixed$Read_Abundance_Sample.t, "gamma")
qqp(soilmixed$Read_Abundance_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#check that mean is above 5 for PQL GLM
mean(soilmixed$Read_Abundance_Sample)
th<-theta.ml(y =soilmixed$Read_Abundance_Sample.t, mu = mean(soilmixed$Read_Abundance_Sample.t, limit = 10))

soilmixed$Read_Abundance_Sample.t <-log(soilmixed$Read_Abundance_Sample.t)
msread <- glm(Read_Abundance_Sample.t ~ SITE + P + NH4 + NO3 + SITE:NH4,  data = soilmixed, family = gaussian(link = "identity"))
stepAIC(msread)
summary(msread)
plot(msread)

########################## Late ###################

# find the right distribution for late roots
rootslate$Read_Abundance_Sample.t <- rootslate$Read_Abundance_Sample + 1
qqp(rootslate$Read_Abundance_Sample.t, "norm")

#LOGNORMAL
qqp(rootslate$Read_Abundance_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(rootslate$Read_Abundance_Sample.t, "Negative Binomial")
qqp(rootslate$Read_Abundance_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(rootslate$Read_Abundance_Sample.t, "Poisson")
qqp(rootslate$Read_Abundance_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(rootslate$Read_Abundance_Sample.t, "gamma")
qqp(rootslate$Read_Abundance_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#check that mean is above 5 for PQL GLM
mean(rootslate$Read_Abundance_Sample)

th<-theta.ml(y =rootslate$Read_Abundance_Sample.t, mu = mean(rootslate$Read_Abundance_Sample.t, limit = 10))

laterich <- glm(Read_Abundance_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3 + TREATMENT*pH 
                + TREATMENT*P + SITE*pH, 
                family =  neg.bin(theta = th), data = rootslate)
stepAIC(laterich)
laterichfinal <- glm(Read_Abundance_Sample.t ~  TREATMENT + SITE + P + NH4 + NO3 + 
                       pH + TREATMENT:SITE + SITE:NH4 + TREATMENT:NO3 + TREATMENT:P + 
                       SITE:pH, 
                     family =  neg.bin(theta = th), data = rootslate)
summary(laterichfinal)
plot(laterichfinal)



# find the right distribution for late soil
soillate$Read_Abundance_Sample.t <- soillate$Read_Abundance_Sample + 1
qqp(soillate$Read_Abundance_Sample.t, "norm")

#LOGNORMAL
qqp(soillate$Read_Abundance_Sample, "lnorm")

#NEG BINOMIAL
nbinom <- fitdistr(soillate$Read_Abundance_Sample.t, "Negative Binomial")
qqp(soillate$Read_Abundance_Sample.t, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])

#POISSON
poisson <- fitdistr(soillate$Read_Abundance_Sample.t, "Poisson")
qqp(soillate$Read_Abundance_Sample.t, "pois", poisson$estimate)

#GAMMA
gamma <- fitdistr(soillate$Read_Abundance_Sample.t, "gamma")
qqp(soillate$Read_Abundance_Sample.t, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#check that mean is above 5 for PQL GLM
mean(soillate$Read_Abundance_Sample)

th<-theta.ml(y =soillate$Read_Abundance_Sample.t, mu = mean(soillate$Read_Abundance_Sample.t, limit = 10))

latelyrich <- glm(Read_Abundance_Sample.t ~ TREATMENT +SITE  + P +NH4 +NO3 + pH + SITE*TREATMENT + SITE*NH4
                  + SITE*P + TREATMENT*NH4 + SITE*NO3 + TREATMENT*NO3 + TREATMENT*pH 
                  + TREATMENT*P + SITE*pH, 
                  family =  gaussian(link = "identity"), data = soillate)
stepAIC(latelyrich)

latelyrichfinal <- glm(formula = Read_Abundance_Sample.t ~ TREATMENT + SITE + P + 
                         NH4 + NO3 + pH + SITE:P + TREATMENT:NH4 + SITE:NO3 + TREATMENT:P, 
                       family = gaussian(link = "identity"), data = soillate)
plot(latelyrichfinal)
summary(latelyrichfinal)
