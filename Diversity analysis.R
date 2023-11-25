#Diversity analysis
#alphadiv####
genus_data <- read.table("genus_data.csv",sep = ",", header = T,row.names=1,stringsAsFactors = FALSE,check.names = FALSE)
genus_data <- genus_data[,colSums(genus_data>1)>10] 
genus_log <- as.data.frame(log10(genus_data+1))

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
env1 <- dplyr::select(genus_metadata, "Disease_Phase","DMM")
env1 <- as.data.frame(unclass(env1), stringsAsFactors = TRUE)
env1$Sample<-rownames(env1)

library(vegan)
Shannon <- diversity(genus_data, index = "shannon")
Simpson <- diversity(genus_data, index = "simpson")
observed_species <- specnumber(genus_data)
Chao1  <- estimateR(genus_data)[2, ]
ACE  <- estimateR(genus_data)[4, ]

alphadata <- data.frame(Shannon,Simpson,observed_species,Chao1,ACE )
alphadata$Sample<-rownames(alphadata)

alphadiv<- alphadata%>%left_join(env1,by = "Sample")

alphadiv <- dplyr::select(alphadiv,"Shannon","Simpson","Disease_Phase", "DMM","Severity")
library(reshape2)
alphadiv.melt <- melt(alphadiv,id.var1 = "Disease_Phase",id.var2 = "DMM",variable.name = "Diversity",value.name = "Value")
alphadiv.melt$Disease_Phase <- factor(alphadiv.melt$Disease_Phase, levels = c("Healthy","Incremental","Critical","Complicated","Convalescent","Long-term follow-up"))
color=c("#4DAF4A", "#FF7F00","#E41A1C", "#377EB8","#984EA3","#A65628")

library(rstatix)
kruskal_test <- alphadiv.melt %>%
  group_by(Diversity) %>%
  kruskal_test(Value ~ Disease_Phase)
dunn_test <- alphadiv.melt %>%
  group_by(Diversity) %>%
  dunn_test(Value ~ Disease_Phase,p.adjust.method = "bonferroni", detailed = FALSE)

library(ggpubr)
mycomparsion <- list(c("Long-term follow-up","Convalescent"),
                     c("Long-term follow-up","Complicated"),
                     c("Long-term follow-up","Critical"),
                     c("Long-term follow-up","Incremental"),
                     c("Healthy", "Incremental"), 
                     c("Healthy","Critical"),
                     c("Healthy","Complicated"),
                     c("Healthy","Convalescent"),
                     c("Healthy","Long-term follow-up"))

F2.1 <- ggboxplot(data=alphadiv.melt, x="Disease_Phase", y="Value", color = "Disease_Phase",
                  palette =  color,bxp.errorbar = T, add = c("jitter","fun"),  short.panel.labs = F)+
  facet_wrap(~ Diversity,scales = "free")+
  stat_compare_means(aes(Disease_Phase), label = "p.format")+
  stat_compare_means(comparisons =mycomparsion, method = "wilcox.test",label = "p.signif",step.increase = 0.05)+
  xlab(NULL)+ylab(NULL)+
  theme1+theme(legend.position ="right",axis.ticks.x = element_blank(),axis.text.x=element_blank())
F2.1
ggsave("alphadiversity-Disease_Phase.pdf",F2.1 ,width = 9,height =6,units = "in",dpi = 300)


alphadiv.melt$DMM <- factor(alphadiv.melt$DMM, levels = c("DMM1","DMM2","DMM3","DMM4","DMM5","DMM6"))
color=c("#FFDC91FF","#20854EFF","#BC3C29FF","#EE4c97FF", "#7876B1FF", "#6F99ADFF")

S2.2 <- ggboxplot(data=alphadiv.melt, x="DMM", y="Value", color = "DMM",
                  palette =  color,bxp.errorbar = T, add = c("jitter","fun"),  short.panel.labs = F)+
  facet_wrap(~ Diversity,scales = "free")+
  stat_compare_means(aes(DMM), label = "p.format")+
  xlab(NULL)+ylab(NULL)+
  theme1+theme(legend.position ="right",axis.ticks.x = element_blank(),axis.text.x=element_blank())
S2.2
ggsave("alphadiversity-DMM.pdf",S2.2 ,width = 9,height =6,units = "in",dpi = 300)

#betadiv####
genus_data <- read.table("genus_data.csv",sep = ",", header = T,row.names=1,stringsAsFactors = FALSE,check.names = FALSE)
genus_data <- genus_data[,colSums(genus_data>1)>10]
genus_log <- as.data.frame(log10(genus_data+1))

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
env2 <- dplyr::select(genus_metadata, "Disease_Phase","DMM","Age","Gender")
env2$Sample<-rownames(env2)
env2 <- as.data.frame(unclass(env2), stringsAsFactors = TRUE)


library(vegan)
bray<-vegdist(genus_log, method="bray")
rownames(env2)<-rownames(as.matrix(bray))

# adjust for age and gender 
library(aPCoA)
result<-aPCoA(bray~ Age + Gender,env2,maincov="Disease_Phase",drawEllipse=TRUE,drawCenter=F,
              pch=19,cex=0.9,lwd=1.2,col=c("#4DAF4A", "#FF7F00","#E41A1C", "#377EB8","#984EA3","#A65628"))

ggsave(result,"aPCoA.pdf",width = 8,height = 8,units = "in",dpi = 300)

aPCoA.score <- data.frame(result$plotMatrix)
sample_point <- data.frame(aPCoA.score[1:2])
names(sample_point)[1:2] <- c('PCoA1', 'PCoA2')
sample_point$Sample=rownames(sample_point)
sample_point<- sample_point%>%left_join(env2,by = "Sample")

library(ggplot2)
sample_point$Disease_Phase <- factor(sample_point$Disease_Phase, levels = c("Healthy","Incremental","Critical","Complicated","Convalescent","Long-term follow-up"))
color=c("#4DAF4A", "#FF7F00","#E41A1C", "#377EB8","#984EA3","#A65628")
F2.2 <- ggplot(sample_point,aes(PCoA1, PCoA2)) +
  geom_point(aes(color=Disease_Phase), size = 2.5) + 
  stat_ellipse(aes(fill = Disease_Phase), geom = 'polygon',
               level = 0.95, alpha = 0.1, show.legend = T) +
  scale_color_manual(values =color) +
  scale_fill_manual(values = color) + 
  labs(x = paste('PCoA1: ',33.57, '%'), y = paste('PCoA2: ',9.35, '%'))+
  theme1+theme(legend.position="right")
F2.2
ggsave("pcoa_Disease_Phase.pdf",F2.2 ,width = 9,height =6,units = "in",dpi = 300)

sample_point$DMM <- factor(sample_point$DMM, levels = c("DMM1","DMM2","DMM3","DMM4","DMM5","DMM6"))
color=c("#FFDC91FF","#20854EFF","#BC3C29FF","#EE4c97FF", "#7876B1FF", "#6F99ADFF")
S2.3 <- ggplot(sample_point,aes(PCoA1, PCoA2)) +
  geom_point(aes(color=DMM), size = 2.5) + 
  stat_ellipse(aes(fill = DMM), geom = 'polygon',
               level = 0.95, alpha = 0.1, show.legend = T) +
  scale_color_manual(values =color) +
  scale_fill_manual(values = color) + 
  labs(x = paste('PCoA1: ',33.57, '%'), y = paste('PCoA2: ',9.35, '%'))+
  theme1+theme(legend.position="right")
S2.3
ggsave("pcoa_DMM.pdf",F2.2 ,width = 9,height =6,units = "in",dpi = 300)

#permutational multivariate analysis of variance (PERMANOVA)
genus.adonis<-adonis2(genus_log ~Disease_Phase+Age+Gender,data=env2, permutations=999,dist = 'bray', by="margin")
genus.adonis

genus.adonis<-adonis2(genus_log ~DMM+Age+Gender,data=env2, permutations=999,dist = 'bray', by="margin")
genus.adonis

#pairwise PERMANOVA
library(pairwiseAdonis)
result<- pairwise.adonis(x=genus_log, factors=env2$Disease_Phase, sim.function = "vegdist",sim.method = "bray",p.adjust.m = "BH",reduce = NULL, perm = 999)
result                        
#write.csv(result,"pairwise.adonis.csv")

#W*d test
source("Wd.R")
genus.Wdtest<- WdS.test(bray, env2$Disease_Phase, nrep = 999) 
#pairwise W*d test
result<-Tw2.posthoc.tests(bray, env2$Disease_Phase, nrep=999)
result 
#write.csv(result,"pairwise.Wdtest.csv")

#permutational analyses of multivariate dispersions (PERMDISP)
dispersion <- betadisper(bray, group=env2$Disease_Phase)

genus.permutest<-permutest(dispersion,permutations=999)
genus.permutest

#pairwise PERMDISP
result <-permutest(dispersion, pairwise = TRUE, permutations = 999)
result 
#write.csv(result,"pairwise.PERMDISP.csv")