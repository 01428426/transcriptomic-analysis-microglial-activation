library(qs)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(cluster)
library(clustree)

results_path <- "/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/"
obj=readRDS('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/combined_object.rds')


DefaultAssay(obj) <- "integrated"

obj <- FindNeighbors(obj, reduction = "umap", dims = 1:2)
obj <- FindClusters(object = obj,
                    resolution = c(0.005,0.075, 0.01,0.015,0.02,0.03, 0.04, 0.05,0.06,0.07,0.08,0.09, 0.1, 0.2))

qs::qsave(obj,"/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/obj_clusters_20pcs.qs")

clust=clustree(obj@meta.data, prefix = "integrated_snn_res.")
ggsave(filename="clustree_new.pdf",
       plot =clust,
       width = 210,
       height = 297,
       units = "mm")

save.image("filename.rda")