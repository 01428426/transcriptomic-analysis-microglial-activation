
library(RColorBrewer)

combined_GO$star=ifelse(combined_GO$`FDR`<=0.01, '**', 
                        ifelse(combined_GO$`FDR`<=0.05,'*',''))
library("data.table")
order=setDT(combined_GO)[ , .N, keyby =Description] 
merge=merge(combined_GO,order, by='Description')

ggplot(merge, aes(as.factor(Cluster),reorder(Description, Ratio), label = star, ylab='genesets')) + 
  geom_tile(aes(fill = Ratio)) + geom_text() +
  scale_fill_gradient2(low = "white", high = "blue")+
  theme(axis.text.x = element_text(angle = 300, vjust = 0, hjust=0))+
  labs(fill = "Enrichment Ratio")+
  ylab("Functional pathways")+
  xlab("Clusters")+
  scale_x_discrete(position="top")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

results1=merge
