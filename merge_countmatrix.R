setwd("/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/2019_Early_embryo/")
directory <- "/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/2019_Early_embryo/"
sampleFiles <- grep(".count",list.files(directory),value=TRUE)
n = length(sampleFiles)     
count.matrix = read.table(file = sampleFiles[1],header=F,sep="\t",row.names = 1)
colnames(count.matrix) = sampleFiles[1]
for (i in 2:n){
  new.data = read.csv(file = sampleFiles[i], header=T, sep="\t",row.names = 1)
  colnames(new.data) = sampleFiles[i]
  count.matrix = merge(count.matrix,new.data,by=0,all=T)
  rownames(count.matrix)=count.matrix$Row.names
}
count.matrix[is.na(count.matrix)]=0
count.matrix<- count.matrix[-(1:5),-(1:(n-1))]
count.matrix<- count.matrix[!(rowSums(count.matrix)==0),]
write.csv(count.matrix,"count.matrix.csv")
