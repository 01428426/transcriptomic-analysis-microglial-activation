library(dplyr)

###### ENRICHMENT FUNCTION
enrichment_custom <- function(genes, reference, genesets, adj = "fdr", verbose = FALSE) {
  tab <- lapply(1:length(genesets), function(i) {
    if (verbose == TRUE) {
      cat("processing term", i, names(genesets)[i], "\n")
    }
    
    ra <- length(genes) / length(reference)
    geneset_size <- length(genesets[[i]])
    expect <- geneset_size * ra
    
    reference <- reference[!reference %in% genes]
    RinSet <- sum(reference %in% genesets[[i]])
    RninSet <- length(reference) - RinSet
    GinSet <- sum(genes %in% genesets[[i]])
    enrichment_ratio <- round(GinSet / expect, 2)
    overlap_genes <- paste(intersect(genes, genesets[[i]]), collapse = ",")
    
    GninSet <- length(genes) - GinSet
    fmat <- matrix(c(GinSet, RinSet, GninSet, RninSet),
                   nrow = 2,
                   ncol = 2, byrow = F
    )
    colnames(fmat) <- c("inSet", "ninSet")
    rownames(fmat) <- c("genes", "reference")
    fish <- fisher.test(fmat, alternative = "greater")
    pval <- as.numeric(format(fish$p.value, format = "e", digits = 2))
    inSet <- RinSet + GinSet
    pct_overlap <- round((GinSet / inSet) * 100, 2)
    res <- c(GinSet, inSet, geneset_size, pct_overlap, enrichment_ratio, pval, overlap_genes)
    res
  })
  rtab <- do.call("rbind", tab)
  rtab <- data.frame(as.vector(names(genesets)), rtab)
  rtab <- rtab[order(rtab[, 4]), ]
  colnames(rtab) <- c(
    "TermID",
    "genes", "all", "geneset_size",
    "pct_overlap", "enrichment_ratio",
    "pval", "overlap_genes"
  )
  tab.out <- rtab %>%
    dplyr::mutate_at(c(
      "genes", "all", "geneset_size",
      "pct_overlap", "enrichment_ratio",
      "pval"
    ), as.character) %>%
    dplyr::mutate_at(c(
      "genes", "all", "geneset_size",
      "pct_overlap", "enrichment_ratio",
      "pval"
    ), as.numeric) %>%
    dplyr::mutate(padj = p.adjust(pval), method = adj) %>%
    dplyr::mutate(padj = as.numeric(format(padj, format = "e", digits = 2)))
  tab.out <- tab.out %>%
    dplyr::select(-overlap_genes, everything(), overlap_genes)
  
  return(tab.out)
}


#####
background_gene_all <- lapply(background_gene_all_l, function(x){x <- x %>%
  dplyr::filter(Cell_proportion > 0.01) %>%
  pull(Markers) %>%
  as.character()})
background_gene_all <- Reduce(union, background_gene_all)

reference=background_gene_all

markers_list <- lapply(markers_list, function(x) {x <- x %>%
  filter(avg_log2FC >0.5 & p_val_adj<=0.05 ) %>%
  pull(gene) %>% as.character()})


# Enrichment top 100
enrich_100=list()
for(i in 1:length(markers_list)){
  enrich_100[[i]]=enrichment_custom(genes=markers_list[[i]],reference,genesets=final_top100)
}

saveRDS(enrich_100,'top_100_enrichment_6.rds')

enrich_100[[1]]
enrich_100[[2]]
enrich_100[[3]]#DAM
enrich_100[[4]] #DAM
enrich_100[[5]]#PVM
enrich_100[[6]]

names(enrich_100)<-c('0','1','2','3','4','5')
library(data.table)

enriched2=as.data.frame.list(enrich_100, nrow=40)
library(data.table)
enriched2= rbindlist(enrich_100)
enriched2$cluster=NA
enriched2$cluster[1:8]='0'
enriched2$cluster[9:16]='1'
enriched2$cluster[17:24]='2'#DAM
enriched2$cluster[25:32]='3'#DAM
enriched2$cluster[33:40]='4'
enriched2$cluster[41:48]='5'

levels(enriched2$TermID)[levels(enriched2$TermID)=='Core micro 3 (Patir et al.)'] <- "Core microglial genes"
levels(enriched2$TermID)[levels(enriched2$TermID)=='Homeostatic (Keren Shaul et al.)'] <- "Homeostatic genes"
levels(enriched2$TermID)[levels(enriched2$TermID)=='aging'] <- "Aging genes"

enriched2$star=ifelse(enriched2$padj >0.05, '', ifelse(
  enriched2$padj <= 0.05 & enriched2$padj > 0.01, '*','**'))

ggplot(enriched2, aes(cluster,reorder(TermID,enrichment_ratio), fill = enrichment_ratio, label = star, ylab='genesets')) + 
  geom_tile() + geom_text() +
  scale_fill_gradient2(low = "white", high = "red")+
  theme(axis.text.x = element_text(angle = 300, vjust = 0, hjust=0))+
  labs(fill = "Enrichment Ratio")+
  ylab("genesets")+
  scale_x_discrete(position="top")

write.csv(enrich_100[[3]],'overlap_DAM1_6.csv')
write.csv(enrich_100[[4]],'overlap_DAM2_6.csv')
write.csv(enrich_100[[5]],'overlap_PVM_6.csv')
write.csv(enriched2,'enrichment_clusters_aim12_6.csv')

