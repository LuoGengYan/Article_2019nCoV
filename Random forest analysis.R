#Random forest analysis
library(randomForest)
library(cowplot)
library(caret)
library(pROC)

genus_RPM <- read.table("genus_RPM.csv",sep = ",", header = T,row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
genus_RPM$Sample<-rownames(genus_RPM)

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
env1 <- dplyr::select(genus_metadata, "Age","Gender","Comorbidity","DMM","Pseudotime")
env1<-subset(env1, env1$Pseudotime %in% c('1','2'))
env1$Sample<-rownames(env1)
env1 <- as.data.frame(unclass(env1), stringsAsFactors = TRUE)

al1_RPM<- env1 %>%left_join(genus_RPM,by = "Sample")
row.names(al1_RPM) <- al1_RPM$Sample

otu_RPM <-al1_RPM[,-c(1:6)]
otu_RPM<- otu_RPM[,colSums(otu_RPM>1)>3]
otu_log <- as.data.frame(log10(otu_RPM+1))
otu_log$Sample<-rownames(otu_log)


otu<- otu_log %>%left_join(env1,by = "Sample")
row.names(otu) <- otu$Sample
otu<-otu[,-c(132)]
for (i in names(otu)[c(133:136)]){otu[,i] <- as.factor(otu[,i])}

#The data is divided into the training set and the verification set by 3:1
set.seed(123)
inTrain <- createDataPartition(y = otu$Pseudotime,p = 0.25,list = F)
validation <- otu[inTrain,]
train <- otu[-inTrain,]
table(train$Pseudotime)
table(validation$Pseudotime)

#Optimize parameters ntree and mtry of RF model
n <- length(names(train))
set.seed(123)
min = 100
num = 0
for (i in 1:(n-1)){
  mtry_fit <- randomForest(Pseudotime~., data = train, mtry = i)
  err <- mean(mtry_fit$err.rate)
  print(err)
  if(err < min) {    
    min = err     
    num = i }
}
print(min)
print(num)
#[1] 127
ntree_fit <- randomForest(Pseudotime~.,data = train,mtry = 127,ntree = 2000)
ntree_fit
pdf("ntree-2 vs 1.pdf",height = 5,width = 5)
plot(ntree_fit)
dev.off()
#write.table(ntree_fit$err.rate,'ntree-2 vs 1.tsv',sep = '\t',quote = F)

#RF model
train.rf <- randomForest(Pseudotime ~ ., data = train, importance = TRUE, proximity = TRUE,mtry = 127,ntree = 2000)
print(train.rf)

#Random Forest Cross-Valdidation for feature selection, the number of features is 29.

#Multiple cross verification
result <- replicate(10, rfcv(train[-ncol(train)], train$Pseudotime, cv.fold = 10, step = 1.5), simplify = FALSE)
result[[1]]$n.var
error.cv <- sapply(result, "[[", "error.cv")
write.table(error.cv,'error.cv-2 vs 1.tsv',sep = '\t',quote = F)

pdf("CV_error-2 vs 1.pdf",height = 5,width = 5)
matplot(result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type = "l",
        lwd = c(2, rep(1, ncol(error.cv))), col = 1, lty = 1, log = "x",
        xlab = "Number of variables", ylab = "Cross-validation error")
dev.off()

#Variable Importance Top 20
pdf("varImpPlot-2 vs 1.pdf",height=5,width=10)
varImpPlot(train.rf, sort = TRUE, n.var = min(20, nrow(train.rf$importance)))
dev.off()

importance_otu <-data.frame(importance(train.rf))
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseGini, decreasing = TRUE), ]
head(importance_otu)
otu_select <- rownames(importance_otu)[1:4]
name_3 <- c(otu_select, "Pseudotime")

train_yc <- train[ ,name_3]
val_yc <- validation[,name_3]

#Cross-validate 10 fold
folds <- createFolds(y = train_yc[,1],k = 10)
fc1 <- as.numeric()
mod_pre1<-as.numeric()
for(i in 1:10){
  fold_test1 <- train_yc[folds[[i]],]
  fold_train1 <- train_yc[-folds[[i]],]
  model1 <- randomForest(Pseudotime ~ .,data = fold_train1,proximity = T,importance = T)
  model_pre1 <- predict(model1,newdata = fold_test1,type="prob")
  fc1 <- append(fc1,as.numeric(fold_test1$Pseudotime))
  mod_pre1 <- append(mod_pre1,model_pre1[,2])
}
df1 <- cbind(fc1,as.numeric(mod_pre1))
write.table(df1,"train_2 vs 1.tsv",sep = '\t',quote = F)

#ROC for train 
pdf("train_2 vs 1.pdf",height = 5,width = 5)
mycol <- c("slateblue","#4DBBD5FF","#000000")
a <- plot.roc(df1[,1],df1[,2],
              smooth = F,
              lwd = 2,
              ylim = c(0,1),
              xlim = c(1,0),
              legacy.axes = T,
              main = "Random Forest",
              ci = TRUE,
              col = mycol[2])
ciobj <- ci.se(a,specificities = seq(0, 1, 0.01))
plot(ciobj, type = "shape", col = "#4DBBD5FF")
legend.name <- c(paste("Critical vs Incremental",sep = " "),paste("AUC","=",round(a$auc,3),sep = " "),
                 paste("95%CI =",round(a$ci[1],3),"-",round(a$ci[3],3),sep = " "))
legend("bottomright", 
       legend = legend.name,
       col = mycol[2:3],
       lwd = 2,
       bty = "n")
dev.off()

#Validation
fc2 <- as.numeric()
mod_pre2 <- as.numeric()
model_FR2 <- randomForest(Pseudotime ~ .,data = train_yc,proximity = T,importance =T )
model_FR_pre2 <- predict(model_FR2,newdata = val_yc,type = "prob")
fc2 <- append(fc2,as.numeric(val_yc$Pseudotime))
mod_pre2 <- append(mod_pre2,model_FR_pre2[,2])
df2 <- cbind(fc2,as.numeric(mod_pre2))
write.table(df2,"validation_2 vs 1.tsv",sep = '\t',quote = F)

#ROC for validation
pdf("validation_2 vs 1.pdf",height = 5,width = 5)
b<- plot.roc(df2[,1],df2[,2],
             smooth = F,
             lwd = 2,
             ylim = c(0,1),
             xlim = c(1,0),
             legacy.axes = T,
             main = "Random Forest",
             ci = TRUE,
             col = mycol[2])
ciobj <- ci.se(b,specificities = seq(0, 1, 0.01))
plot(ciobj, type="shape", col = "#4DBBD5FF")
legend.name <- c(paste("Critical vs Incremental",sep = " "),paste("AUC","=",round(b$auc,3),sep = " "),
                 paste("95%CI=",round(b$ci[1],3),"-",round(b$ci[3],3),sep = " "))
legend("bottomright", 
       legend = legend.name,
       col = mycol[2:3],
       lwd = 2,
       bty = "n")
dev.off()
