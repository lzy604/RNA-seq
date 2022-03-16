#install.packages
#install.packages("ggpubr")

library("ggpubr")

#caculate the correlation,pearson mean one batch,spearman mean different batch
cor(x, y, method = c("pearson", "kendall", "spearman"))

#plot
ggscatter(my_data, x = "mpg", y = "wt", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")
#my_data mean your data
#x,y mean which your plot from the data
ggscatter(data, x, y, combine = FALSE, merge = FALSE,
          color = "black", fill = "lightgray", palette = NULL, shape = 19,
          size = 2, point = TRUE, rug = FALSE, title = NULL, xlab = NULL,
          ylab = NULL, facet.by = NULL, panel.labs = NULL,
          short.panel.labs = TRUE, add = c("none", "reg.line", "loess"),
          add.params = list(), conf.int = FALSE, conf.int.level = 0.95,
          fullrange = FALSE, ellipse = FALSE, ellipse.level = 0.95,
          ellipse.type = "norm", ellipse.alpha = 0.1,
          ellipse.border.remove = FALSE, mean.point = FALSE,
          mean.point.size = ifelse(is.numeric(size), 2 * size, size),
          star.plot = FALSE, star.plot.lty = 1, star.plot.lwd = NULL,
          label = NULL, font.label = c(12, "plain"), font.family = "",
          label.select = NULL, repel = FALSE, label.rectangle = FALSE,
          cor.coef = FALSE, cor.coeff.args = list(), cor.method = "pearson",
          cor.coef.coord = c(NULL, NULL), cor.coef.size = 4, ggp = NULL,
          show.legend.text = NA, ggtheme = theme_pubr(), ...)
#corrplot
cor_matrix <- cor (as.matrix(fpkm.all))
corrplot(corr=cor_matrix, add = TRUE, type = "lower", method = "number",
         col =brewer.pal(10,"RdGy"))
png("cor.fpkm.png",width = 8,height = 8)
corrplot.mixed(cor_matrix,tl.pos="lt",cl.lim = c(0,1),number.font=0.02,col)
