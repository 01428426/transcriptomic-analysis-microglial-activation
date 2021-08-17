library(cluster)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(MAST)
library(qs)

obj_umap=qread('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/obj_clusters_20pcs.qs')


DefaultAssay(obj_umap) <- "integrated"
Idents(object = obj_umap) <- "integrated_snn_res.0.015"

markers=readRDS('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/markers.rds')


library(dplyr)


top5 <- markers%>% 
  group_by(cluster) %>% top_n(5, avg_log2FC)


vector=top5$gene
print(vector)

# set.seed(111)
# random.cells <- sample(x = colnames(obj_umap@assays$integrated), size = 10000, replace = F)
# print(random.cells)
# 
h <- obj_umap@meta.data %>%
  group_by(integrated_snn_res.0.015)%>%
  group_by(diagnosis)%>%
  sample_frac(size = 0.2, replace=F)

random.cells=as.character(h$barcode)
length(random.cells)
random.cells

FeaturePlot(obj_umap, features = c("TREM2"), pt.size=0.01)
DoHeatmap(obj_umap, cells=random.cells,features = c("TREM2","CD58","ERAP2","GNLY","S100A12","TNF","IL1","CD68","TYROBP"), size=2, group.by= c("integrated_snn_res.0.015", "diagnosis"))  
