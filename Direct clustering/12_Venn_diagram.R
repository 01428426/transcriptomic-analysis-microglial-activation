
install.packages('VennDiagram')
suppressPackageStartupMessages(library(VennDiagram))


 ###
 
 
 ### Five Sets
 
 de1=read.table(file='CTR_vs_diagnosisAD1.tsv',sep = '\t', header = TRUE, fill = TRUE)
 
 de2=read.table(file='CTR_vs_diagnosisAD2.tsv',sep = '\t', header = TRUE, fill = TRUE)
 
 de3=read.table(file='CTR_vs_diagnosisAD3.tsv',sep = '\t', header = TRUE, fill = TRUE)
 
 de4=read.table(file='CTR_vs_diagnosisAD4.tsv',sep = '\t', header = TRUE, fill = TRUE)
 
 de5=read.table(file='CTR_vs_diagnosisAD5.tsv',sep = '\t', header = TRUE, fill = TRUE)
 
 de6=read.table(file='CTR_vs_diagnosisAD6.tsv',sep = '\t', header = TRUE, fill = TRUE)
 
 
 bonf = list(`cluster 0`=gene1,`cluster 1`=gene2,`cluster 2`=gene3,`cluster 3`=gene4,
             `cluster 4`=gene5,`cluster 5`=gene6)
 
 ### Six Sets
 library(venn)
 
 dataForVennDiagram <-bonf
 
 #vennDiagramColors <- c('#EA4335', '#FBBC05', '#34A853', '#4285F4', 'orchid3')
 vennDiagramColors <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3", 'cyan') # TODO: cyan is not good enough
 vennFilename <- 'Venn_Diagram_6Sets.png'
 
 venn(dataForVennDiagram, zcolor = "style", opacity = 0.25, cexil = 1, ilcs=1, sncs=1,cexsn = 0.1, ellipse = FALSE, box=FALSE)
 
 library(ggVennDiagram)
 venn.plot <-ggVennDiagram(bonf,label_alpha = 0,stroke_color = "white", digits=0,label = 'count',label_color = "white",show_intersect=FALSE)
 
 
 
 
 