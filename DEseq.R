#DESeq2差异分析
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("tximport")
library(DESeq2)
library(stringr)
#Create a data frame with columns as specified below containing the necessary information for the design matrix (sampleName, fileName, condition).
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/20200707_epsc/count/")
directory <- "/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/20200707_epsc/count/"
sampleFiles <- grep(".count",list.files(directory),value=TRUE)
sampleCondition <- sub(".count","",sampleFiles)
sampleCondition <- c(rep("TdEPS",4),rep("EPSC2C",4),rep("EPSC8C",4),rep("ESC",2))
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
#import the data and format the data for analysis
directory <- "/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/20200707_epsc/count/"
sampleFiles <- grep(".count",list.files(directory),value=TRUE)
n = length(sampleFiles)     
count.matrix = read.table(file = sampleFiles[1],header=F,sep="\t")
colnames(count.matrix) = c("gene",sampleFiles[1])
rownames(count.matrix)=count.matrix$gene
for (i in 2:n){
  new.data = read.csv(file = sampleFiles[i], header=T, sep="\t")
  colnames(new.data) = c("gene",sampleFiles[i])
  rownames(new.data)=new.data$gene
  count.matrix = merge(count.matrix,new.data,all=T)
}
count.matrix[is.na(count.matrix)]=0
rownames(count.matrix)=count.matrix$gene
count.matrix<- count.matrix[-(1:5),-1]
#deseq2
TdEPSvsEPS2C <- DESeqDataSetFromMatrix(countData = count.matrix[,1:8],colData=sampleTable[1:8,],design=~condition)
TdEPSvsEPS2C <- TdEPSvsEPS2C[ rowSums(counts(TdEPSvsEPS2C)) > 1, ]   #过滤low count数据
nrow(TdEPSvsEPS2C)
TdEPSvsEPS2C <- DESeq(TdEPSvsEPS2C)
TdEPSvsEPS2C.results <- results(TdEPSvsEPS2C)
summary(TdEPSvsEPS2C.results)
write.csv(TdEPSvsEPS2C.results,file="TdEPSvsEPS2C.results.csv")

TdEPSvsEPS2C <- DESeqDataSetFromMatrix(countData = count.matrix[,1:8],colData=sampleTable[1:8,],design=~condition)
TdEPSvsEPS2C <- TdEPSvsEPS2C[ rowSums(counts(TdEPSvsEPS2C)) > 1, ]   #过滤low count数据
nrow(TdEPSvsEPS2C)
TdEPSvsEPS2C <- DESeq(TdEPSvsEPS2C)
TdEPSvsEPS2C.results <- results(TdEPSvsEPS2C)
summary(TdEPSvsEPS2C.results)
write.csv(TdEPSvsEPS2C.results,file="TdEPSvsEPS2C.results.csv")
