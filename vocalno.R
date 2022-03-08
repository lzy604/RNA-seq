vocanp_plot=function(dd,genes=NULL,filename,cutoffpvalue,cutofflogfc,dir){
  dd$Gene.symbol=rownames(dd)
  p=ggplot(dd, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = group_rnaseq)) +
    scale_color_manual(values = c("blue", "grey","red")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  if (length(genes)>0){
    p=p+geom_label_repel(data =dd[genes,],aes(label = Gene.symbol,
                                              fill = group_rnaseq), color = 'black',
                         box.padding = 1, label.padding = 0.24,
                         arrow = arrow(length = unit(0.02, "npc")),
                         size = 3.5) +
      scale_fill_manual(values = setNames(c("lightblue", "lightsalmon",'gray','aliceblue'), levels(dd$group_crispr)))+
      theme(legend.position = "bottom")
  }
  p=p+geom_vline(xintercept = cutofflogfc,linetype="dashed")+
    geom_vline(xintercept = -cutofflogfc,linetype="dashed")+
    geom_hline(yintercept = -log10(cutoffpvalue),linetype="dashed")
  p=p+theme (axis.text.x=element_text(angle = 0,hjust=1,size = 16,face = "bold"),
             axis.text.y=element_text(size = 16,face = "bold"),
             title=element_text(size = 14,face = "bold"),
             axis.title.x=element_text(size = 15,face = "bold"),
             legend.text=element_text(size = 12,face = "bold"),
             legend.title=element_text(size = 12,face = "bold"),
             legend.position="none")
  print(p)
  ggsave(filename = paste0(dir,'/',filename,'.pdf'),width = 6,height = 6)
}