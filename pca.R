#pca
# https://www.r-bloggers.com/computing-and-visualizing-pca-in-r/
# https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

#package
#install.packages("devtools", repo="http://cran.us.r-project.org")
library(devtools)
#install.packages("ggplot2")
library(ggplot2)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(sva)

#read file
eps.fpkm <- read.csv("/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/20200707_epsc/fpkm/eps.fpkm.csv",header = T,row.names = 1)
embryo.fpkm <- read.table("/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/2019_Early_embryo/gene.FPKM.combined.tab",header = T)
embryo.fpkm<- embryo.fpkm[!duplicated(embryo.fpkm$GeneId),]
row.names(embryo.fpkm) <- embryo.fpkm$GeneId
embryo.fpkm <- embryo.fpkm[,13:15]
fpkm.all <- merge(eps.fpkm,embryo.fpkm,by = 0,all = T)
row.names(fpkm.all) <- fpkm.all$Row.names
fpkm.all <- fpkm.all[,-1]
fpkm.all[is.na(fpkm.all)] <-0
fpkm.all <- fpkm.all[which(rowSums(fpkm.all)>0),]
write.csv(fpkm.all,"/Users/liziyi/Documents/gaolab.data/EPS/NGS/RNA-seq/20200707_epsc/fpkm/fpkm.all.csv")

#PCA
log2fpkm.all <- log2(fpkm.all+1)#log transformed
pca.log2fpkm.all <- prcomp(t(log2fpkm.all), center = TRUE, scale. = TRUE)#pca,remember if you use the sample to do the pca,you need to transpose
#plot(pca.log2fpkm.all, type = "l")
#summary(pca.log2fpkm.all)
fpkm_pc_df <- as.data.frame(pca.log2fpkm.all$x)
fpkm_pc_df$status <- c(rep("cellline",14),rep("embryo",3))
fpkm_pc_df$Condition <- c(rep("TdEPS",4),rep("EPSC2C",4),rep("EPSC8C",4),rep("ESC",2),"2cell","4cell","8cell")
fpkm_pc_df$Condition <- c(rep("TdEPS",4),rep("EPSC2C",4),rep("EPSC8C",4),rep("ESC",2),"2cell","4cell","8cell")

#pca_PLOT
pdf(file = 'log2fpkm_pca.pdf')
ggplot(fpkm_pc_df, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)
dev.off()

#batch remove and pca
batch.info <- data.frame(1:17,row.names = colnames(fpkm.all))
colnames(batch.info)  <- "samples"
batch.info$batch <- c(rep(1,14),rep(2,3))
batch.info$type <- c(rep("cellline",14),rep("embryo",3))
log2fpkm.all <- log2fpkm.all[which(rowSums(log2fpkm.all) >0),]
combat_exp <- ComBat(dat = as.matrix(log2fpkm.all), batch = batch.info$batch,  mod=NULL)

pca.combat_exp.log2fpkm.all <- prcomp(t(combat_exp), center = TRUE, scale. = TRUE)
combat.log2fpkm.all.pc <- as.data.frame(pca.combat_exp.log2fpkm.all$x)
combat.log2fpkm.all.pc$status <- fpkm_pc_df$status
combat.log2fpkm.all.pc$Condition  <- fpkm_pc_df$Condition

pdf(file = 'log2fpkm_combat_pca.pdf')
ggplot(combat.log2fpkm.all.pc, aes(x = PC1, y = PC2, color = Condition)) +
  geom_point(size=4)
dev.off()

