#DMM clusters####
library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(ggplot2)
##Read raw RPM file
genus_RPM <- read.table("genus_RPM.csv",sep = ",", header = T,row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
genus_RPM <- genus_RPM[,colSums(genus_RPM>1)>10]

count<-as.matrix(genus_RPM)

#the fitting effect
fit<-lapply(1:10,dmn,count=count,verbose=TRUE)
lplc<-sapply(fit,laplace)
aic<-sapply(fit,AIC)
bic<-sapply(fit,BIC)
plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
+lines(aic, type="b", lty = 2)
+lines(bic, type="b", lty = 3)

S2.1 <-plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
ggsave("alphadiversity-DMM.pdf",S2.1 ,width = 9,height =6,units = "in",dpi = 300)



best<-fit[[which.min(unlist(lplc))]]
mixturewt(best)
heatmapdmn(count, fit[[1]], best, 30)
DMM <- apply(mixture(best), 1, which.max)
write.csv(DMM, file="DMM_6clusters.csv")

for(k in seq(ncol(fitted(best)))){
  d<-melt(fitted(best))
  colnames(d)<-c("OTU","cluster","value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    filter(abs(value) > quantile(abs(value), 0.8))  
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() + 
    labs(title = paste("Top drivers: community type", k))
  print(p)
}

#Kaplan-Meier curves
library(autoReg)
library(dplyr)
library(survival)
library(survminer) 
genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
mytable <- dplyr::select(genus_metadata,"Gender","Age","Comorbidity","Days_post_symptom_onset","Pseudotime","DMM","Clinical_outcome")
mytable <- subset(mytable , mytable $Clinical_outcome %in% c('Severe','Deceased'))

fit <- survfit( Surv(Days_post_symptom_onset,Clinical_outcome) ~ DMM,  data = mytable)
summary(fit)
ggsurvplot(fit,  pval = TRUE,surv.median.line = "hv",conf.int = F,risk.table = TRUE,
           legend.title = "DMM",
           legend.labs = c("DMM1", "DMM2","DMM3", "DMM4","DMM5", "DMM6"),
           palette =  c("#20854EFF","#FFDC91FF","#BC3C29FF","#EE4c97FF", "#7876B1FF", "#6F99ADFF"),
           xlim = c(0,66), 
           xlab = "Days post symptom onset")

#log-rank test
survdiff(Surv(Days_post_symptom_onset,Clinical_outcome) ~ DMM, data = mytable)
#pairwise test
pairwise_survdiff(Surv(Days_post_symptom_onset,Clinical_outcome) ~ DMM,data = mytable,p.adjust.method = "BH")        


#adjusted for age, gender and comorbidities
psModel<-glm(Clinical_outcome~Gender+Age+Comorbidity,family=binomial(link="logit"),data=mytable) 
mytable$ps=predict(psModel,type="response")

mytable$IPTW<-ifelse(mytable$DMM=="DMM1",1/mytable$ps,1/(1-mytable$ps))

fit.IPTW<- survfit(Surv(Days_post_symptom_onset,Clinical_outcome) ~ DMM, 
                   weights=mytable$IPTW,
                   data = mytable)
summary(fit.IPTW)

p<- ggsurvplot(fit.IPTW,  pval = TRUE,surv.median.line = "hv",conf.int = F,risk.table = F,
           legend.title = "DMM",
           legend.labs = c("DMM1", "DMM2","DMM3", "DMM4","DMM5", "DMM6"),
           palette =  c("#20854EFF","#FFDC91FF","#BC3C29FF","#EE4c97FF", "#7876B1FF", "#6F99ADFF"),
           xlim = c(0,66), 
           xlab = "Days post symptom onset")
p
ggsave("DMM-KM.pdf",p,width = 8,height =6,units = "in",dpi = 300)
