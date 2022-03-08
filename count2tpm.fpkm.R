#read file
setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/20200707_epsc/count/")
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
write.csv(count.matrix,"count.matrix.csv")
#read gene length
gene.length <- read.csv("mm10_gene_len.txt",header=T,sep="\t")
row.names(gene.length)<- gene.length$gene
#TPM
count.with.length <- merge(count.matrix,gene.length,by=0)
rpk <- 1000*count.with.length[,(2:15)]/count.with.length$length
denom <- (colSums(rpk))/1e6
count2tpm <- data.frame(row.names=count.with.length$Row.names)
for (i in 1:ncol(rpk) ){
  tpm1 = as.data.frame(rpk[,i]/denom[i],row.names=count.with.length$Row.names)
  count2tpm <- cbind(count2tpm,tpm1)
}
colnames(count2tpm) <- colnames(count.matrix)
write.csv(count2tpm,"all.count2tpm.csv",sep = "\t")
#FPKM
count.with.length <- merge(count.matrix,gene.length,by=0)
sum.count <- (colSums(count.with.length[,2:15]))/1e6
rpm <-  data.frame(row.names=count.with.length$Row.names)
for (i in 2:(ncol(count.with.length)-3) ){
  tpm1 = as.data.frame(count.with.length[,i]/sum.count[i-1],row.names=count.with.length$Row.names)
  rpm <- cbind(rpm,tpm1)
}
colnames(rpm) <- colnames(count.matrix)
count2rpkm <- 1000*rpm/count.with.length$length
