library(qs)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)
library(dplyr)


sce_human=qread('/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/gerrits_smith/processed_sce/Micro_sce.qs')

### HUMAN DATASET 


#### filter human dataset based on mouse cutoff
sce_subset <- sce_human[ , sce_human$total_counts >=400]

#switch the rownames from ensembl id from gene name 
#by the following commands and then convert it to seurat object.
sce_subset<- sce_subset[!is.na(SummarizedExperiment::rowData(sce_subset)$gene), ]
sce_subset <- sce_subset[!duplicated(SummarizedExperiment::rowData(sce_subset)$gene), ]
rownames(sce_subset) <- SummarizedExperiment::rowData(sce_subset)$gene

sce_subset$species <- 'human'


seu<- CreateSeuratObject(counts = SingleCellExperiment::counts(sce_subset), 
                   meta.data = as.data.frame(SingleCellExperiment::colData(sce_subset), project = "human"))

seu[["RNA"]] <- AddMetaData(seu[["RNA"]], as.data.frame(SummarizedExperiment::rowData(sce_subset)),col.name = NULL)


### CLUSTERS CORRESPONDING TO OPTIMIZED RESOLUTION
obj=qread("/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/round_2/obj_clusters_20pcs.qs")
obj <- obj[ , obj$total_counts >=400]

seu@meta.data$clusters=obj@meta.data$integrated_snn_res.0.08

####merge the two seurat objects
mouse_seu=qs::qread("/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/round_2/mouse_seu.qs")
  
df1 = data.frame(name=colnames(seu@assays$RNA))
df2=data.frame(name=colnames(mouse_seu@assays$RNA))
names=merge(df1,df2, by='name')
seu=seu[which(rownames(seu@assays$RNA)%in% names[,1]),]
mouse_seu=mouse_seu[which(rownames(mouse_seu@assays$RNA)%in% names[,1]),]

qsave(mouse_seu,'/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/ML/mouse_seu_filtered.qs')
sce.seurat<- merge(mouse_seu, seu, add.cell.ids = c("human", "mouse"), project = "human_mouse")
rownames(sce.seurat@assays$RNA)
gene_list=list(rownames(sce.seurat@assays$RNA))
lenght(gene_list)
saveRDS(gene_list,'/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/ML/gene_liste.rds'))

### INTEGRATION

# split the dataset into a list of  seurat objects (by species/manifest)
seu.list <- SplitObject(sce.seurat, split.by = "manifest")
# normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 500)
})
# select features that are repeatedly variable across datasets for integration
# run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = seu.list)
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

#perform integration
anchors <- FindIntegrationAnchors(object.list = seu.list, anchor.features = features, reduction = "rpca")
combined_human_mouse <- IntegrateData(anchorset =anchors, dims=1:20, k.weight=50)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(combined_human_mouse) <- "integrated"

# Run the standard workflow for visualization and clustering
combined_human_mouse <- ScaleData(combined_human_mouse, verbose = FALSE)

qs::qsave(combined_human_mouse,"/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/round_2/combined_mouse_human_manifest.qs")

