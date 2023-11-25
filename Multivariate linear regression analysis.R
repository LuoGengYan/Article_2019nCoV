#Multivariate linear regression test
library(autoReg)
library(dplyr)
library(survival)
library(survminer) 

Clini <- read.table("clinical-indicators.csv",sep = ",", header = T,row.names=1,stringsAsFactors = FALSE,check.names = FALSE)
Clini <- as.data.frame(log2(Clini+1))
Clini$Sample<-rownames(Clini)

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
genus_metadata$Sample<-row.names(genus_metadata)
env1 <- dplyr::select(genus_metadata,"Sample","PatientID","Age","Gender","Comorbidity","Pseudotime","DMM","Clinical_outcome")
env1 <- subset(env1 , env1 $Pseudotime %in% c("1", "2", "3", "4"))                    

env1$Mortality[env1$Clinical_outcome=="Deceased"]=1 
env1$Mortality[env1$Clinical_outcome!="Deceased"]=0 
env1$Mortality=as.factor(env1$Mortality)

lm_t =Clini%>%left_join(env1,by = "Sample")%>%subset(PatientID != "na")


lm_t_res=data.frame(beta=0,p=0)
for(i in 1:ncol(lm_t[,c(1:21)])){
  fm=as.formula(paste("lm_t[,i]~",paste(names(lm_t)[c(30)],collapse = "+"),sep = ""))
  tmp=lm(fm,data = lm_t2) %>% summary()
  lm_t_res[i,]=c(tmp$coefficients[2,1],tmp$coefficients[2,4])
}
rownames(lm_t_res)=names(lm_t)[c(1:21)]
lm_t_res$p.adj=p.adjust(lm_t_res$p,method = "fdr",n = nrow(lm_t_res))

write.csv(lm_t_res,"lm_t_res-Mortality.csv")


lm_t1 <- subset(lm_t,lm_t$"Pseudotime" == "1")
lm_t1_res=data.frame(beta=0,p=0)
for(i in 1:ncol(lm_t1[,c(1:21)])){
  fm=as.formula(paste("lm_t1[,i]~",paste(names(lm_t1)[c(30)],collapse = "+"),sep = ""))
  tmp=lm(fm,data = lm_t2) %>% summary()
  lm_t1_res[i,]=c(tmp$coefficients[2,1],tmp$coefficients[2,4])
}
rownames(lm_t1_res)=names(lm_t1)[c(1:21)]
lm_t1_res$p.adj=p.adjust(lm_t1_res$p,method = "fdr",n = nrow(lm_t1_res))

write.csv(lm_t1_res,"lm_t1_res-Mortality.csv")

lm_t2 <- subset(lm_t,lm_t$"Pseudotime" == "1")
lm_t2_res=data.frame(beta=0,p=0)
for(i in 1:ncol(lm_t2[,c(1:21)])){
  fm=as.formula(paste("lm_t2[,i]~",paste(names(lm_t2)[c(30)],collapse = "+"),sep = ""))
  tmp=lm(fm,data = lm_t2) %>% summary()
  lm_t2_res[i,]=c(tmp$coefficients[2,1],tmp$coefficients[2,4])
}
rownames(lm_t2_res)=names(lm_t2)[c(1:21)]
lm_t2_res$p.adj=p.adjust(lm_t2_res$p,method = "fdr",n = nrow(lm_t2_res))

write.csv(lm_t2_res,"lm_t2_res-Mortality.csv")


##correlation
Clini <- read.table("clinical-indicators.csv",sep = ",", header = T,row.names=1,stringsAsFactors = FALSE,check.names = FALSE)
Clini <- as.data.frame(log2(Clini+1))
Clini$Sample<-rownames(Clini)

genus_RPM <- read.table("genus_RPM.csv",sep = ",", header = T,row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
genus_RPM <- genus_RPM[,colSums(genus_RPM>1)>10]
genus_RPM$Sample<-rownames(genus_RPM)

genus_Clini <- Clini%>%left_join(genus_RPM ,by = "Sample")%>%subset(Acinetobacter.baumannii != "na")

vars=c("WBC","HGB","PLT","NEUT","LYMPH","NLR","DD","ALT","AST", "LDH","TG","CRP","SAA",
       "PCT","IL_6", "T_lymphocyte", "CD4_T_lymphocyte", "CD8_T_lymphocyte", "CD4_CD8", "B_cell", "NK_cell")    

df_tmp=genus_Clini[,vars]
genus <- genus_Clini[,tmp1] %>% t()

genus <- log10(genus+1)

pval.df=data.frame(matrix(0,length(vars),nrow(genus)))
r.df=data.frame(matrix(0,length(vars),nrow(genus)))
names(pval.df)=rownames(genus)
names(r.df)=rownames(genus)
rownames(pval.df)=vars
rownames(r.df)=vars

####correlation
for(i in vars){
  for (j in rownames(genus)){
    tmp=cor.test(df_tmp[,i] %>% as.numeric(),
                 genus[j,] %>% as.numeric())
    pval.df[i,j]=tmp$p.value
    r.df[i,j]=tmp$estimate
  }
}

p.adjust.foo= function(x) p.adjust(x,method = "fdr")
pval.df = apply(pval.df,1,p.adjust.foo) %>% t()

nrow(pval.df)
sig.thres = 0.01

pval.df.sig = pval.df[rowSums(pval.df<sig.thres)>0,colSums(pval.df<sig.thres)>0]
r.df.sig =r.df[na.omit(rownames(pval.df.sig)),na.omit(colnames(pval.df.sig))]
pval.df.sig =pval.df.sig[na.omit(rownames(r.df.sig)),na.omit(colnames(r.df.sig))]

pre_genus=0
for (i in 1:ncol(pval.df.sig)) {
  tmp=which(pval.df.sig[,i]<0.01)
  if(max(abs(r.df.sig[tmp,i]))>0.3){
    pre_genus[i]=colnames(pval.df.sig)[i]
  }
}
pre_genus=na.omit(pre_genus)
pre_genus=pre_genus[pre_genus !=0]

pval.df.sig=pval.df.sig[,pre_genus]
r.df.sig=r.df.sig[,pre_genus]

colnames(r.df.sig)
if("WBC" %in% rownames(r.df.sig)){
  sig=colnames(r.df.sig)[order(r.df.sig["WBC",],decreasing = T)]
}else{sig=colnames(r.df.sig)}


bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01)) 
r.df.sig=t(r.df.sig)
pval.df.sig=t(pval.df.sig)

pheatmap(r.df.sig[sig,],
         color = colorRampPalette(c("#13457E", "white", "#A71E31"))(length(bk)),
         #angle_col = "45",
         border=F,
         display_numbers = matrix(ifelse((pval.df.sig[sig,]<=sig.thres & pval.df.sig[sig,]>0.01),"", 
                                         ifelse((pval.df.sig[sig,]<=0.01 & pval.df.sig[sig,]>0.001),"**",
                                                ifelse(pval.df.sig[sig,]<=0.001,"***",""))), nrow(pval.df.sig)),
         cluster_rows = T,
         cluster_columns = T,
         fontsize = 10,
         cellwidth = 15, cellheight =10,
         cluster_cols = T,breaks = bk,border_color = "grey60",
         na_col = "grey90"
)
save_pheatmap_pdf <- function(x, Genus.Clini, width=12, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(Genus.Clini))
  pdf(Genus.Clini, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
