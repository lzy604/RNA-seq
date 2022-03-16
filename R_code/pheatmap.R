source("http://biocoundctor.org/biocLite.R")
biocLite("pheatmap")

library(pheatmap)

pdf("all_heatmap_cluster.pdf")
pheatmap(tmp, #表达数据
         show_rownames = F,# 显示行名
         cluster_cols = F,#对列是否做cluster
         cluster_rows = F,
         show_colnames = T,# 显示列名
         scale = "row", #对行标准化
         color =colorRampPalette(c("#8854d0", "#ffffff","#fa8231"))(100), # 热图基准颜色
         border_color = 'NA')
dev.off()



# color = colorRampPalette(c("navy", "white", "firebrick3"))(100) 蓝红

# clustering_distance_rows = "correlation"，表示行聚类使用皮尔森相关系数聚类，当然也可以自定义如drows = dist(test, method = "minkowski")；clustering_distance_rows = drows
#cluster_row = FALSE#表示行不聚类
#legend = FALSE#表示右侧图例不显示
#display_numbers = TRUE#表示在热图中格子显示对应的数字，在那种横纵轴数目比较小是时候可用，比如样本间相关系数聚类
#number_format = "\%.1e"#当显示数字时数字的显示方式
#cellwidth = 15, cellheight = 12#表示热图中小方格的宽度和高度
#fontsize = 8#表示热图中字体显示的大小
#filename = "test.pdf"#表示直接就保存成test.pdf图片了
#labels_row#可以自己定义横轴的显示字符，默认上图是基因名
#main#类似title啦
#gaps_col#产生一个间隔，就像有些文章中的那种分类后每个分类都有一个间隔
#annotation_col = data.frame(ClassGene = factor(paste0('Cluster',cutree(clust$tree_col,10)))，
