library(qs)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)


sce=qread('/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/gerrits_smith/processed_sce/Micro_sce.qs')

sce <- sce[!is.na(SummarizedExperiment::rowData(sce)$gene), ]
sce <- sce[!duplicated(SummarizedExperiment::rowData(sce)$gene), ]
rownames(sce) <- SummarizedExperiment::rowData(sce)$gene

seu <- CreateSeuratObject(counts = SingleCellExperiment::counts(sce),
                          meta.data = as.data.frame(SingleCellExperiment::colData(sce)))

seu[["RNA"]] <- AddMetaData(seu[["RNA"]],
                            as.data.frame(SummarizedExperiment::rowData(sce)),
                            col.name = NULL)

saveRDS(seu,'/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/round_2/seu.rds')
