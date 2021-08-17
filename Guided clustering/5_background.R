library(cluster)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(qs)
####
# updated 1/31/2020 to accommodate V3.1
# updated 2/4/2020 to output "NA" for genes not detected in certain subgroups
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(parallel::mclapply(genes, mc.cores = 1, calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    results = parallel::mclapply(list, mc.cores = 1, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    #combined = do.call("rbind", results)
    return(results)
  }
}
calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}
####

seu_micro_pc=qread('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/clara_aim12/obj_clusters_manifest_20pcs.qs')
DefaultAssay(seu_micro_pc) <- "integrated"
Idents(object = seu_micro_pc) <- "integrated_snn_res.0.005"


system.time(background_gene_all_l <- PrctCellExpringGene(object = seu_micro_pc,
                                                         genes = rownames(seu_micro_pc),
                                                         group.by = "integrated_snn_res.0.005"))
#it took about 25 min
background_gene_all <- lapply(background_gene_all_l, function(x){x <- x %>%
  dplyr::filter(Cell_proportion > 0.25) %>%
  pull(Markers) %>%
  as.character()})
background_gene_all <- Reduce(union, background_gene_all)
save(background_gene_all_l, background_gene_all,
     file = "/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/clara_aim12/background_genes_exprsd_0.25prct_cell.rda")