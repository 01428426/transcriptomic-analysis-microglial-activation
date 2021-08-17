library(qs)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)
library(SingleCellExperiment)
# 
sce=qread('/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/gerrits_smith/processed_sce/Micro_sce.qs')
seu=qread('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/obj_clusters_20pcs.qs')

idx <- which(seu$barcode %in% sce$barcode)
idx <- seu$barcode[idx]
seu_subset <- subset(seu, cells = names(idx))
sce_subset <- sce[ , idx]
all(sce_subset$barcode == seu_subset$barcode)

sce_subset$cluster=seu_subset@meta.data$integrated_snn_res.0.015


sce_1=sce_subset[ , which(sce_subset$cluster=='0')]

sce_2=sce_subset[ , which(sce_subset$cluster=='1')]

sce_3=sce_subset[ , which(sce_subset$cluster=='2')]

sce_4=sce_subset[ , which(sce_subset$cluster=='3')]

sce_5=sce_subset[ , which(sce_subset$cluster=='4')]

sce_6=sce_subset[ , which(sce_subset$cluster=='5')]



qsave(sce_1,"sce_1.qs")
qsave(sce_2,"sce_2.qs")
qsave(sce_3,"sce_3.qs")
qsave(sce_4,"sce_4.qs")
qsave(sce_5,"sce_5.qs")
qsave(sce_6,"sce_6.qs")
