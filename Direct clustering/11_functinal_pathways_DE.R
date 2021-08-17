library(qs)


library(dplyr)
library(sctransform)
library(ggplot2)


de1=read.table(file='CTR_vs_diagnosisAD1.tsv',sep = '\t', header = TRUE, fill = TRUE)

de2=read.table(file='CTR_vs_diagnosisAD2.tsv',sep = '\t', header = TRUE, fill = TRUE)

de3=read.table(file='CTR_vs_diagnosisAD3.tsv',sep = '\t', header = TRUE, fill = TRUE)

de4=read.table(file='CTR_vs_diagnosisAD4.tsv',sep = '\t', header = TRUE, fill = TRUE)

de5=read.table(file='CTR_vs_diagnosisAD5.tsv',sep = '\t', header = TRUE, fill = TRUE)

de6=read.table(file='CTR_vs_diagnosisAD6.tsv',sep = '\t', header = TRUE, fill = TRUE)

de1_filter = de1 %>%
  filter(logFC>0.5)%>%
  filter(padj<0.05)

gene1=as.character(de1_filter$gene)


de2_filter = de2 %>%
  filter(logFC>0.5)%>%
  filter(padj<0.05)

gene2=as.character(de2_filter$gene)



de3_filter = de3 %>%
  filter(logFC>0.5)%>%
  filter(padj<0.05)

gene3=as.character(de3_filter$gene)



de4_filter = de4 %>%
  filter(logFC>0.5)%>%
  filter(padj<0.05)

gene4=as.character(de4_filter$gene)



de5_filter = de5 %>%
  filter(logFC>0.5)%>%
  filter(padj<0.05)

gene5=as.character(de5_filter$gene)

de6_filter = de6 %>%
  filter(logFC>0.5)%>%
  filter(padj<0.05)

gene6=as.character(de6_filter$gene)

write.table(gene1,col.names = F, row.names = F , quote = FALSE,"gene1.txt")
write.table(gene2,col.names = F, row.names = F , quote = FALSE,"gene2.txt")
write.table(gene3,col.names = F, row.names = F , quote = FALSE,"gene3.txt")
write.table(gene4,col.names = F, row.names = F , quote = FALSE,"gene4.txt")
write.table(gene5,col.names = F, row.names = F , quote = FALSE,"gene5.txt")
write.table(gene6,col.names = F, row.names = F , quote = FALSE,"gene6.txt")
write.table(reference,col.names = F, row.names = F , quote = FALSE,"reference.txt")

