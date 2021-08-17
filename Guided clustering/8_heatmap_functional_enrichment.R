combined_GO <- read_delim("combined_GO.csv", 
                          +     ";", escape_double = FALSE, trim_ws = TRUE)

library(RColorBrewer)

combined_GO$star=ifelse(combined_GO$`P Value`<0.05 & combined_GO$`P Value` >= 0.01, '*', '**')
library("data.table")
order=setDT(combined_GO)[ , .N, keyby =Description] 
merge=merge(combined_GO,order, by='Description')
merge$cluster=as.factor(merge$cluster)

ggplot(merge, aes(cluster,reorder(Description, N), label = star)) + 
  geom_tile(aes(fill = Ratio)) + geom_text() +
  scale_fill_gradient2(low = "white", high = "blue")+
  labs(fill = "Enrichment Ratio")+
  ylab("Functional pathways")+
  xlab("Clusters")+
  scale_x_discrete(position="top")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

results2=merge
