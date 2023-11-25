
#Pairwise differential abundance analysis####
#####DESeq2####
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
##Read raw count file
genus_count <- read.table("genus_count.csv",sep = ",",row.names = 1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
genus_count <- genus_count[ ,colSums(genus_count) >= 1.5*ncol(genus_count)]
genus_count$Sample <-  row.names(genus_count)

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
env1 <- dplyr::select(genus_metadata, "PatientID","Age","Gender","Pseudotime")
env1$Sample<-rownames(env1)
env1 <- as.data.frame(unclass(env1), stringsAsFactors = TRUE)


count_data <- genus_count%>%left_join(env1,by = "Sample") 
row.names(count_data)<- count_data$Sample
genus_count <- as.data.frame(t(count_data[,1:141]+1))

#Read reference gene info file
mart_export <- read.table("Tax.csv",sep = ",", header = T,stringsAsFactors = FALSE,check.names = FALSE)
unique_mart <- subset(mart_export, duplicated(mart_export$Tax_ID) == FALSE)
rownames(unique_mart) <- unique_mart$Tax_name


env1 <- count_data[,-c(1:141)]
env1$Sample <- row.names(env1) 
env1$Pseudotime<-factor(env1$Pseudotime)

pseudotime <- c("1", "2", "3", "4", "5")
folder <- "DESeq2_results-COVID-vs-Control"
if (!dir.exists(folder)){
  dir.create(folder)
}

#Compare each pseudotime with control samples
for (i in 1:length(pseudotime)) {
  output_folder <- file.path(folder, paste("Control_vs_COVID", pseudotime[i], sep = ''))
  if (!dir.exists(output_folder)){
    dir.create(output_folder)
  }
  
  
  col_data <- subset(env1, env1$Pseudotime %in% c('0', pseudotime[i]))
  rownames(col_data) <- col_data$Sample
  col_data$pseudotime <- as.factor(as.character(col_data$Pseudotime))
  

  selected_genus_count <- genus_count[, as.character(col_data$Sample)]
  
  
  dds_counts <- DESeqDataSetFromMatrix(countData =round(selected_genus_count) , colData = col_data, design = ~ Pseudotime)
  dds_counts <- estimateSizeFactors(dds_counts)
  
  
  dds <- DESeq(dds_counts, betaPrior = FALSE)
  res <- results(dds, independentFiltering = TRUE, alpha = 0.05)
  res_sorted <- res[order(res$padj), ]
  
  
  res_sorted$tax <- unique_mart[rownames(res_sorted), ]$Tax_name
  write.table(res_sorted, file = file.path(output_folder, paste("DESeq2result_taxnames_Control_vs_COVID_", pseudotime[i], ".txt", sep = '')), 
              sep="\t", quote=FALSE)
  top_gene <- rownames(res_sorted)[1]
} 

#Code to generate volcano plot###

plot_data <- data.frame()


for (i in 1:length(pseudotime)) {
  input_folder <- file.path(folder, paste("Control_vs_COVID", pseudotime[i], sep = ''))
  result <- read.csv(file.path(input_folder, paste("DESeq2result_taxnames_Control_vs_COVID_", pseudotime[i], ".txt", sep = '')), 
                     sep = '\t')
  result <- subset(result, result$padj < 0.05 & result$baseMean > 100 & abs(result$log2FoldChange) > 0.5)
  result$pseudotime <- pseudotime[i]
  result <- rownames_to_column(result, "Tax_name")
  plot_data <- rbind(plot_data, result)
  
}
write.csv(plot_data,file = "DESeq2result-0.csv")

plot_data$direction <- ifelse(plot_data$log2FoldChange < 0, "down", "up")


Tax_frequency <- table(plot_data$Tax_name)
plot_data$frequency <- as.vector(Tax_frequency[as.character(plot_data$Tax_name)])


tax_to_label <- c("Morococcus", "Aggregatibacter", "Fusobacterium",
                  "Campylobacter", "Alloprevotella","Neisseria", "Dialister", "Tannerella",
                  "Streptococcus", "Lachnoanaerobaculum", "Haemophilus", "Capnocytophaga",
                  "Granulicatella","Rothia", "Porphyromonas","Parvimonas","Prevotella","Veillonella")
tax_to_label2 <- c("Acinetobacter","Actinomyces", "Aspergillus","Betacoronavirus", "Candida","Enterococcus", "Corynebacterium", "Herbaspirillum",  "Mycoplasma",
                   "Citrobacter","Escherichia","Klebsiella","Simplexvirus","Pseudomonas","Ralstonia","Cupriavidus","Staphylococcus","Dolosigranulum")


plot_data1 <- subset(plot_data, plot_data$padj < 0.05 & plot_data$baseMean > 100 & abs(plot_data$log2FoldChange) >= 3)
#Code to plot 
cairo_ps(file = "DESeq2result-0.ps", width = 8, height =10)
p <- ggplot(data = plot_data1, aes(x=pseudotime, y=log2FoldChange, fill = direction))
p <- p + geom_point(position=position_jitter(width=0.1), aes(size=-1*log10(padj), alpha=frequency), pch=21, color="black", stroke=0.1)
p <- p + geom_text_repel(data = subset(plot_data1, plot_data1$tax %in% tax_to_label & abs(plot_data1$log2FoldChange) > 1.5), aes(x=pseudotime, y=log2FoldChange, label=tax), 
                         color="black", size=3.5, nudge_x = -0.5)
p <- p + geom_text_repel(data = subset(plot_data1, plot_data1$tax %in% tax_to_label2& abs(plot_data1$log2FoldChange) > 1.5), aes(x=pseudotime, y=log2FoldChange, label=tax), 
                         color="black", size=3.5, nudge_x = -0.5)
p <- p + scale_color_manual(values=c('#004C99', '#990000')) + scale_fill_manual(values=c('#004C99', '#990000'))
p <- p + xlab("") + ylab("Log fold change") + scale_x_discrete(labels = c("Incremental vs Healthy", "Critical vs Healthy", "Complicated vs Healthy", "Convalescent vs Healthy", "Long-term follow-up vs Healthy"))
p <- p + theme_bw() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=12), 
                            plot.title = element_text(size=20, face="bold", hjust = 0.5), axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()

#####Maaslin2####
library(Maaslin2)
genus_RPM <- read.table("genus_RPM.csv",sep = ",", header = T,row.names = 1,stringsAsFactors = FALSE, check.names = FALSE)
genus_RPM <- genus_RPM[,colSums(genus_RPM>1)>10]
genus_RPM <-genus_RPM+1
genus_RPM$Sample<-rownames(genus_RPM)

genus_metadata <- read.table("genus_metadata.csv",sep = ",",row.names=1, header = T,stringsAsFactors = FALSE,check.names = FALSE)
env1 <- dplyr::select(genus_metadata, "PatientID","Age","Gender","Pseudotime")
env1$Sample<-rownames(env1)
env1 <- as.data.frame(unclass(env1), stringsAsFactors = TRUE)

al1_RPM<- genus_RPM %>%left_join(env1,by = "Sample")
row.names(al1_RPM) <- al1_RPM$Sample

#compare with Healthy control (Pseudotime 0)
fit_data <-Maaslin2(
  input_data = al1_RPM[,1:125], 
  input_metadata = al1_RPM[,-c(1:125)], 
  normalization = "NONE",transform = "LOG",analysis_method = "LM", max_significance = 0.05,
  output = "Maaslin2-0", plot_heatmap = TRUE, plot_scatter = F,heatmap_first_n = 50,
  fixed_effects = c("Pseudotime"),
  random_effects = c("Age","Gender"),
  reference = "Pseudotime,0")

#compare amomg samples from COVID-19 patients
al1_RPM<-subset(al1_RPM, al1_RPM$Pseudotime %in% c('1','2','3','4'))
fit_data <-Maaslin2(
  input_data = al1_RPM[,1:125], 
  input_metadata = al1_RPM[,-c(1:125)], 
  normalization = "NONE",transform = "LOG",analysis_method = "LM", max_significance = 0.05,
  output = "Maaslin2-1", plot_heatmap = TRUE, plot_scatter = F,heatmap_first_n = 50,
  fixed_effects = c("Pseudotime"),
  random_effects = c("Age","Gender","PatientID"),
  reference = "Pseudotime,1")

