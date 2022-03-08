library(VennDiagram)
venn.diagram(x=list("diff.gene"=diff.gene,"diff.tr"=diff.tr),imagetype='png',
             filename = "diff.gene.tr.png" ,fill=c("dark red","green"),
             main="difference between genes and transcripts")