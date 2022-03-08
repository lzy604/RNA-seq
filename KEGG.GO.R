#https://zhuanlan.zhihu.com/p/35510434
#https://www.jianshu.com/p/7361c729c79e

# bioconductor包安装
source("https://bioconductor.org/biocLite.R")
biocLite("DOSE") 
biocLite("topGO")
biocLite("clusterProfiler")
biocLite("pathview")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("org.Hs.eg.db")

##加载需要的R包
library(DOSE)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

library(topGO)
library(clusterProfiler)
library(pathview)

##数据载入及转化
##需要将输入的基因格式转为clusterProfiler可分析的格式，通过功能函数bitr实现各种ID之间的转化。通过keytypes函数可查看所有支持及可转化类型
keytypes(org.Hs.eg.db) 
keytypes(org.Mm.eg.db)

#基因list文件读入
data <- read.table("gene",header=FALSE) #单列基因名文件
data$V1 <- as.character(data$V1) #需要character格式，然后进行ID转化
#将SYMBOL格式转为ENSEMBL和ENTERZID格式 
test1 = bitr(data$V1, fromType="SYMBOL", toType=c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
head(test1,2)

#OUT
#SYMBOL         ENSEMBL ENTREZID
#1    AASDH ENSG00000157426   132949
#2   ABCB11 ENSG00000073734     8647

#GO分析
##3.1 groupGO 富集分析
ggo <- groupGO(gene = test1$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)

##3.2 enrichGO 富集分析
ego_ALL <- enrichGO(gene = test1$ENTREZID, 
                    universe = names(geneList), #背景基因集
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 0.01, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 0.05,
                    readable = TRUE) #Gene ID 转成gene Symbol ，易读
head(ego_ALL,2)
#ONTOLOGY         ID                                                Description
#GO:0002887       BP GO:0002887 negative regulation of myeloid leukocyte mediated immunity
#GO:0033004       BP GO:0033004                negative regulation of mast cell activation
#GeneRatio  BgRatio      pvalue  p.adjust    qvalue                   geneID Count
#GO:0002887     2/121 10/11461 0.004706555 0.7796682 0.7796682              CD300A/CD84     2
#GO:0033004     2/121 10/11461 0.004706555 0.7796682 0.7796682              CD300A/CD84     2
#其中：ONTOLOGY：CC  BP  MF 
#GO ID: Gene Ontology数据库中唯一的标号信息
#Description ：Gene Ontology功能的描述信息
#GeneRatio：差异基因中与该Term相关的基因数与整个差异基因总数的比值
#BgRation：所有（ bg）基因中与该Term相关的基因数与所有（ bg）基因的比值
#pvalue: 富集分析统计学显著水平，一般情况下， P-value < 0.05 该功能为富集项
#p.adjust 矫正后的P-Value
#qvalue：对p值进行统计学检验的q值
#geneID：与该Term相关的基因
#Count：与该Term相关的基因数

#3.4 输出结果及图像
##3.4.1 输出结果
write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)

##3.4.2 绘制图形
##可视化--点图
dotplot(ego_MF,title="EnrichmentGO_MF_dot")#点图，按富集的数从大到小的
##可视化--条形图
barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF")#条状图，按p从小到大排，绘制前20个Term

##四、KEGG分析
##4.1 候选基因进行通路分析
kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 1)
head(kk,2)
#ID                                      Description GeneRatio  BgRatio
#hsa04750 hsa04750 Inflammatory mediator regulation of TRP channels      5/53  97/7387
#hsa04020 hsa04020                        Calcium signaling pathway      6/53 182/7387
#pvalue   p.adjust    qvalue                              geneID Count
#hsa04750 0.0006135305 0.08589427 0.0807277             40/3556/3708/5608/79054     5
#hsa04020 0.0018078040 0.12654628 0.1189345        493/1129/2066/3707/3708/4842     6

#4.2 富集结果及图形展示
##4.2.1 结果输出文件
write.csv(summary(kk),"KEGG-enrich.csv",row.names =FALSE)

##4.2.1 KEGG 气泡图
dotplot(kk,title="Enrichment KEGG_dot")
