#Longitudinal GAMM analysis 
# Dynamics of opportunistic species in COVID-19 patients (Cohort 1)
pacman::p_load(tidyverse, dplyr, reshape2,ggsci,ggplot2)

species_RPM <- read.table("species_RPM.csv",sep = ",", header = T,stringsAsFactors = FALSE, row.names = 1,check.names = FALSE)
species_RPM <-species_RPM[,colSums(species_RPM>1)>3]
species_RPM$Sample <-row.names(species_RPM)

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
env1 <- dplyr::select(genus_metadata,"Cohort","PatientID","Days_post_symptom_onset","Pseudotime")
env1 <- subset(env1, env1$Cohort %in% c("COVID_19"))
env1 <-env1[,-1]%>% dplyr::rename(day="Days_post_symptom_onset")
env1$Sample<-rownames(env1)
env1 <- as.data.frame(unclass(env1), stringsAsFactors = TRUE) 

al1_species<-env1%>%left_join(species_RPM,by = "Sample")
row.names(al1_species)<-al1_species$Sample
al1_species<-al1_species[,-3]

mydata <- melt(al1_species,id.vars=c("PatientID","day","Pseudotime"),variable.name = "pathogen",value.name = "RPM")
mydata <- subset(mydata,mydata$pathogen%in% c("SARS-CoV-2",	"HSV-1","Candida albicans",	"Candida glabrata",
                                              "Acinetobacter baumannii","Corynebacterium striatum",
                                              "Enterococcus faecalis","Enterococcus faecium",
                                              "Herbaspirillum huttiense","Klebsiella pneumoniae",
                                              "Staphylococcus epidermidis","Streptococcus","Veillonella"))
mydata$pathogen=factor(mydata$pathogen,levels=c("SARS-CoV-2",	"HSV-1","Candida albicans",	"Candida glabrata",
                                                "Acinetobacter baumannii","Corynebacterium striatum",
                                                "Enterococcus faecalis","Enterococcus faecium",
                                                "Herbaspirillum huttiense","Klebsiella pneumoniae",
                                                "Staphylococcus epidermidis","Streptococcus","Veillonella"))
                          
colour=c( "#FF0000FF", "#FF00CCFF", "#FF9900FF", "#FFCC00FF", 
          "#0000CCFF", "#666600FF",   
          "#6699FFFF", "#99CCFFFF", 
          "#358000FF","#00FFFFFF",   
          "#9900CCFF","#FFCCCCFF", "#999999FF")
library(ggbump)
mydata$day[mydata$day > 66]='70'
mydata$day<-as.numeric(mydata$day)
S.6 <- 
  ggplot(data=mydata,aes(x=day,y=log10(RPM+1),colour=pathogen))+
  facet_wrap(~PatientID,scales = "free",ncol = 5)+
  geom_bump()+ 
  geom_point(size=3)+ 
  scale_colour_manual(values=colour)+
  scale_x_continuous(breaks = integer_breaks()) +
  scale_y_continuous(limit = c(0,7.5))+
  theme+labs(x="Days post symptom onset")
S.6
ggsave("Patient-species.pdf",S.6,width = 17,height =22,units = "in",dpi = 300,limitsize = FALSE)
dev.off()

#longitudinal GAMM analysis (pseudotime 1-4)
library(mgcv) 
library(gamm4)
library(lme4) 
p<-ggplot(mydata, aes(x = Pseudotime, y = mydata$RPM, colour =mydata$pathogen)) + 
  stat_smooth(method="gam", formula = y~s(x,k=4), se=F) + 
  scale_colour_manual(values=colour)+ 
  scale_x_continuous(breaks = seq(1,4,1)) + 
  labs( y = "Log10(RPM+1)",title = "Pseudotime") +xlab(NULL)+
  theme+theme(legend.position ="none",strip.text = element_blank())
p
ggsave("species-gam.pdf",width = 6,height =5,units = "in",dpi = 300,limitsize = FALSE)
dev.off()

dt <- subset(mydata,mydata$pathogen == "SARS-CoV-2")
fml <- "RPM~s(Pseudotime,bs='cr',k=4)"
fit<-gamm(formula(fml),random=list(PatientID=~1),data=dt,family=gaussian)
summary(fit$gam)   
summary(fit$lme)   
pvalue<-summary(fit$gam)[8]
pvalue
anova(fit$gam)

#longitudinal GAMM analysis (different Clinical outcome)
species_RPM <- read.table("species_RPM.csv",sep = ",", header = T,stringsAsFactors = FALSE, row.names = 1,check.names = FALSE)
species_RPM <-species_RPM[,colSums(species_RPM>1)>3]
species_RPM$Sample <-row.names(species_RPM)

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
env1 <- dplyr::select(genus_metadata,"Cohort","PatientID","Days_post_symptom_onset","Clinical_outcome")
env1 <- subset(env1, env1$Cohort %in% c("COVID_19"))
env1 <-env1[,-1]%>% dplyr::rename(day="Days_post_symptom_onset")
env1$Sample<-rownames(env1)
env1 <- as.data.frame(unclass(env1), stringsAsFactors = TRUE) 

al1_species<-env1%>%left_join(species_RPM,by = "Sample")
row.names(al1_species)<-al1_species$Sample
al1_species<-al1_species[,-3]

mydata <- melt(al1_species,id.vars=c("day","Clinical_outcome"),variable.name = "pathogen",value.name = "RPM")
mydata$"Clinical_outcome"=factor(mydata$"Clinical_outcome")

p<-ggplot(mydata, aes(x = day, y = log10(RPM +1), color = Clinical_outcome)) + 
  facet_wrap(~pathogen,ncol =5)+
  stat_smooth(method="gam", formula = y~s(x), se= TRUE) + 
  scale_color_npg() + 
  scale_x_continuous(breaks = seq(0,70,10)) + 
  labs(x = "Days post symptom onset", y = "log10(RPM +1)") + 
  theme1+theme(legend.position ="top")
p

dt <- subset(mydata,mydata$pathogen == "SARS_CoV_2")
fm2 <- RPM ~ Clinical_outcome+s(day, bs = "cr")
model <- gamm4(formula(fm2) , random = ~(1 | PatientID), data =dt)
summary(model$gam)
anova(model$gam)

