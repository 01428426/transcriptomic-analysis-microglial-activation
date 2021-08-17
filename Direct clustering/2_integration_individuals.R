library(Seurat)
#https://satijalab.org/seurat/archive/v3.1/pbmc3k_tutorial.html

sce.seurat=readRDS('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/round_2/seu.rds')
sce.seurat <- sce.seurat[ , sce.seurat$total_counts >=400]
results_path='/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/'

# split the dataset into a list of  seurat objects (by individual) 
seu.list <- SplitObject(sce.seurat, split.by = "individual")
# normalize and identify variable features for each dataset independently
seu.list <- lapply(X = seu.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
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
combined <- IntegrateData(anchorset =anchors, dims=1:20)

# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"

# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)

####
#qr' and 'y' must have the same number of rows error occurs when there are NA values
#solved by removing the cells
####
na_cells <- colSums(as.matrix(is.na(combined@assays$integrated@data)))
to_remove <- names(na_cells)[na_cells>0]
if(length(to_remove) > 0){
  print("Found some cells with NA values. They will be removed. These are the cells:")
  print(to_remove)
  combined@meta.data$keep <- TRUE
  combined@meta.data[to_remove, "keep"] <- FALSE
  combined <- subset(combined, keep)
  combined@meta.data <- combined@meta.data[, -c(which(colnames(combined@meta.data) == "keep"))]
}
####

#Dimensionality reduction
combined <- RunPCA(combined, npcs = 20, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)

#save RDS 
saveRDS(combined, '/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/combined_object.rds')

