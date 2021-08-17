library(dplyr)
library(sctransform)
library(ggplot2)
library(ggrepel)

de1=read.table(file='CTR_vs_diagnosisAD1.tsv',sep = '\t', header = TRUE, fill = TRUE)

de2=read.table(file='CTR_vs_diagnosisAD2.tsv',sep = '\t', header = TRUE, fill = TRUE)

de3=read.table(file='CTR_vs_diagnosisAD3.tsv',sep = '\t', header = TRUE, fill = TRUE)

de4=read.table(file='CTR_vs_diagnosisAD4.tsv',sep = '\t', header = TRUE, fill = TRUE)

de5=read.table(file='CTR_vs_diagnosisAD5.tsv',sep = '\t', header = TRUE, fill = TRUE)

de6=read.table(file='CTR_vs_diagnosisAD6.tsv',sep = '\t', header = TRUE, fill = TRUE)


# Convert directly in the aes()
p <- ggplot(data=de1, aes(x=logFC, y=-log10(padj))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de1, aes(x=logFC, y=-log10(padj))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de1$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de1$diffexpressed[de1$logFC > 0.5 & de1$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de1$diffexpressed[de1$logFC < -0.5 & de1$padj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de1, aes(x=logFC, y=-log10(padj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de1$delabel <- NA
de1$delabel[de1$diffexpressed != "NO"] <- as.character(de1$gene[de1$diffexpressed != "NO"])

plot1= ggplot(data=de1, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()+
  ggtitle('cluster 0')


#############


######
# Convert directly in the aes()
p <- ggplot(data=de2, aes(x=logFC, y=-log10(padj))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de2, aes(x=logFC, y=-log10(padj))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de2$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de2$diffexpressed[de2$logFC > 0.5 & de2$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de2$diffexpressed[de2$logFC < -0.5 & de2$padj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de2, aes(x=logFC, y=-log10(padj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de2$delabel <- NA
de2$delabel[de2$diffexpressed != "NO"] <- as.character(de2$gene[de2$diffexpressed != "NO"])

plot2= ggplot(data=de2, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()+
  ggtitle('cluster 1')




#######

# Convert directly in the aes()
p <- ggplot(data=de3, aes(x=logFC, y=-log10(padj))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de3, aes(x=logFC, y=-log10(padj))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de3$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de3$diffexpressed[de3$logFC > 0.5 & de3$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de3$diffexpressed[de3$logFC < -0.5 & de3$padj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de3, aes(x=logFC, y=-log10(padj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de3$delabel <- NA
de3$delabel[de3$diffexpressed != "NO"] <- as.character(de3$gene[de3$diffexpressed != "NO"])

plot3= ggplot(data=de3, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()+
  xlim(-5,5)+
  ggtitle('cluster 2')


#####
# Convert directly in the aes()
p <- ggplot(data=de4, aes(x=logFC, y=-log10(padj))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de4, aes(x=logFC, y=-log10(padj))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de4$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de4$diffexpressed[de4$logFC > 0.5 & de4$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de4$diffexpressed[de4$logFC < -0.5 & de4$padj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de4, aes(x=logFC, y=-log10(padj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de4$delabel <- NA
de4$delabel[de4$diffexpressed != "NO"] <- as.character(de4$gene[de4$diffexpressed != "NO"])

plot4= ggplot(data=de4, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()+
  ggtitle('cluster 3')

#####
# Convert directly in the aes()
p <- ggplot(data=de5, aes(x=logFC, y=-log10(padj))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de5, aes(x=logFC, y=-log10(padj))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de5$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de5$diffexpressed[de5$logFC > 0.5 & de5$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de5$diffexpressed[de5$logFC < -0.5 & de5$padj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de5, aes(x=logFC, y=-log10(padj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.005), col="red")
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de5$delabel <- NA
de5$delabel[de5$diffexpressed != "NO"] <- as.character(de5$gene[de5$diffexpressed != "NO"])

plot5= ggplot(data=de5, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()+
  ggtitle('cluster 4')



####
# Convert directly in the aes()
p <- ggplot(data=de6, aes(x=logFC, y=-log10(padj))) + geom_point()
# Add more simple "theme"
p <- ggplot(data=de6, aes(x=logFC, y=-log10(padj))) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
de6$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de6$diffexpressed[de6$logFC > 0.5 & de6$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de6$diffexpressed[de6$logFC < -0.5 & de6$padj < 0.05] <- "DOWN"

# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=de6, aes(x=logFC, y=-log10(padj), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
de6$delabel <- NA
de6$delabel[de6$diffexpressed != "NO"] <- as.character(de6$gene[de6$diffexpressed != "NO"])

plot6= ggplot(data=de6, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel()+
  ggtitle('cluster 5')


######
de1$cluster=0
de2$cluster=1
de3$cluster=2
de4$cluster=3
de5$cluster=4
de6$cluster=5

differential=rbind(de1,de2,de3,de4,de5,de6)

plot= ggplot(data=differential, aes(x=logFC, y=-log10(padj), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps=20)+
  facet_grid(vars(cluster), scales = "free")+
  ggtitle('Volcano plots by cluster')


saveRDS(differential,'combined_de.rds')
saveRDS(de1,'de1.rds')
saveRDS(de2,'de2.rds')
saveRDS(de3,'de3.rds')
saveRDS(de4,'de4.rds')
saveRDS(de5,'de5.rds')
saveRDS(de6,'de6.rds')
