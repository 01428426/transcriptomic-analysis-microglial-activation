threshold2=list()
for(col in colnames(threshold)){
  threshold2[[col]]=threshold[,col]
}


markers_l <- readRDS("markers_list.rds")
markers_list <- lapply(markers_list, function(x) {x <- x %>%
  filter(avg_log2FC >= 0.5) %>%
  pull(gene) %>% as.character()})

saveRDS(markers_list,'markers_l2.rds')

genesets <- readRDS("genesets.rds")
genesets_l <- lapply(threshold2, function(dt){dt <- dt[1][[1]]})



external_marker_enrichment_res <- lapply(markers_l, function(genes){
  enrichment_custom(genes = genes, 
                    reference = background_gene_all, 
                    genesets = genesets_l, adj = "fdr", verbose = FALSE)
})

saveRDS(genesets_l,'genesets_l_threshold.rds')
saveRDS(markers_l,'markers_l.rds')

