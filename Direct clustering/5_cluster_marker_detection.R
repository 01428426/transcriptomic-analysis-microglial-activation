
library(cluster)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(MAST)
library(qs)


obj_umap=qread("/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/obj_clusters_20pcs.qs")

Idents(object = obj_umap) <- "integrated_snn_res.0.015"

markers <- FindAllMarkers(obj_umap, 
                            min.pct = 0.25, 
                            logfc.threshold = 0.25,
                            test.use = "wilcox")

head(markers)

marker_list <- split(markers, markers$cluster)

saveRDS(markers,'/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/markers.rds')
saveRDS(marker_list,'/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/markers_list.rds')

  
  