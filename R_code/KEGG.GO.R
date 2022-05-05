#https://zhuanlan.zhihu.com/p/35510434
#https://www.jianshu.com/p/7361c729c79e

# bioconductor install
source("https://bioconductor.org/biocLite.R")
biocLite("DOSE") 
biocLite("topGO")
biocLite("clusterProfiler")
biocLite("pathview")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")

##library
library(DOSE)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

library(topGO)
library(clusterProfiler)
library(pathview)

##data loading
##trasnfer to clusterProfiler file
keytypes(org.Hs.eg.db) 
keytypes(org.Mm.eg.db)

#gene list file
data <- read.table("gene",header=FALSE) #gene name
data$V1 <- as.character(data$V1) #gene id
#SYMBOLto ENSEMBL or ENTERZID
test1 = bitr(data$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
head(test1,2)

#OUT
#SYMBOL         ENSEMBL ENTREZID
#1    AASDH ENSG00000157426   132949
#2   ABCB11 ENSG00000073734     8647

#GO分析
##3.1 groupGO
ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)

##3.2 enrichGO 
ego_ALL <- enrichGO(gene = test1$ENTREZID, 
                    universe = names(geneList), #background
                    OrgDb = org.Hs.eg.db, #organism="human"，OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #or CC  BP  MF
                    pAdjustMethod = "BH", #" holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 0.01, #P
                    qvalueCutoff = 0.05,
                    readable = TRUE) #Gene ID to gene Symbol
head(ego_ALL,2)
#ONTOLOGY         ID                                                Description
#GO:0002887       BP GO:0002887 negative regulation of myeloid leukocyte mediated immunity
#GO:0033004       BP GO:0033004                negative regulation of mast cell activation
#GeneRatio  BgRatio      pvalue  p.adjust    qvalue                   geneID Count
#GO:0002887     2/121 10/11461 0.004706555 0.7796682 0.7796682              CD300A/CD84     2
#GO:0033004     2/121 10/11461 0.004706555 0.7796682 0.7796682              CD300A/CD84     2


#3.4 result and plot
##3.4.1 write.csv
write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)

##3.4.2 plot
##dot
dotplot(ego_MF,title="EnrichmentGO_MF_dot")
##bar
barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF")#top 20 Term

##KEGG
##4.1 
kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'hsa', #organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
#ID                                      Description GeneRatio  BgRatio
#hsa04750 hsa04750 Inflammatory mediator regulation of TRP channels      5/53  97/7387
#hsa04020 hsa04020                        Calcium signaling pathway      6/53 182/7387
#pvalue   p.adjust    qvalue                              geneID Count
#hsa04750 0.0006135305 0.08589427 0.0807277             40/3556/3708/5608/79054     5
#hsa04020 0.0018078040 0.12654628 0.1189345        493/1129/2066/3707/3708/4842     6

#4.2 result and plot
##4.2.1 write.csv
write.csv(summary(kk),"KEGG-enrich.csv",row.names =FALSE)

##4.2.1 KEGG plot
dotplot(kk,title="Enrichment KEGG_dot")
