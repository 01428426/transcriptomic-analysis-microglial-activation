library(scPred)
library(Seurat)
library(SeuratObject)
library(magrittr)
library(dplyr)
library(qs)
library(caret)
library(spatstat)
library(mda)
library(DMwR)
library(pbapply)
library(methods)





######
.make_names <- function(x){
  x <- gsub("\\+", "_plus", x)
  x <- gsub("\\-", "_minus", x)
  x <- make.names(x)
}


.intersectMat <- function(ref, new){
  
  # Get gene names from loadings reference
  refRows <- rownames(ref)
  # Get gene names from new dataset
  newRows <- rownames(new)
  
  if(!(all(newRows %in% refRows) & all(refRows %in% newRows))){ # Subset genes if necesary
    newSub <- new[newRows %in% refRows, ] 
    refSub <- ref[refRows %in% newRows, ]
    newSub <- newSub[match(rownames(refSub), rownames(newSub)), ]
  }else if(!all(newRows == refRows)){ # Order new data according to loadings matrix
    newSub <- newData[match(refRows, newRows), ]
    refSub <- ref
  }else{ # Use data directly if genes match and are ordered
    newSub <- new
    refSub <- ref
  }
  
  return(list(ref = refSub, new = newSub))  
  
}


subsetMatrix <- function(x, s, by.col = TRUE, drop = FALSE, verbose = FALSE, ...){
  
  if(by.col){
    ids <- colnames(x, ...)
    if(is.null(ids) & verbose) stop("No colnames were found!")
    i <- ids %in% s
    if(!any(i) & verbose){message("No matches were found")}
    x <- x[,i, drop = drop]
  }else{
    ids <- rownames(x, ...)
    if(is.null(ids) & verbose) stop("No rownames were found!")
    i <- ids %in% s
    if(!any(i) & verbose){message("No matches were found")}
    x <- x[i, , drop = drop]
  }
  
  x
  
}


getPalette <- function(n){
  
  if(n < 6){
    c("#29BF12", "#00A5CF", "#DE1A1A", "#574AE2", "#FFBF00")
  }else if(n < 9){
    c("#558aa6", "#B1740F", "#D5006A", "#08585A", "#FFFD98", "#9449d2", "#BBBE64", "#D7263D")
  }else if(n < 13){
    c("#943CB4", "#194D44", "yellow", "#5B6DC8", "#3CA437", "#6B244C", "#6ACDC5", "#DE1A1A", "#BBB53E", "#2A297A", "#995533", "#D590DA")
  }else{
    stop("Too many classes")
  }
  
}

###
trainModel_new <- function(object,
                           model = "svmRadial",
                           preProcess = c("center", "scale"),
                           resampleMethod = "cv",
                           number = 5,
                           seed = 66,
                           tuneLength = 10,
                           metric = c("ROC", "PR", "Accuracy", "Kappa"),
                           returnData = FALSE,
                           savePredictions = "final",
                           allowParallel = FALSE,
                           reclassify = NULL
){
  
  
  # Validations -------------------------------------------------------------
  
  # Check class
  if(!is(object, "Seurat") | is(object, "scPred")){
    stop("object must be 'Seurat' or 'scPred'")
  }
  
  
  if(is(object, "Seurat")){
    seurat_object <- object
    object <- get_scpred(object)
    
    if(is.null(object))
      stop("No features have been determined. Use 'getFeatureSpace()' function")
    
    object_class <- "Seurat"
    
  }else{
    object_class <- "scPred"
  }
  
  
  if(is.null(reclassify)){
    classes <- names(object@features)
  }else{
    classes <- reclassify
  }
  metric <- match.arg(metric)
  reduction <- object@reduction
  
  # Train a prediction model for each class
  cat(crayon::green(cli::symbol$record, " Training models for each cell type...\n"))
  
  
  
  if(length(classes) == 1){
    modelsRes <-  .trainModel(classes[1],
                              object,
                              model,
                              reduction,
                              preProcess,
                              resampleMethod,
                              tuneLength,
                              seed,
                              metric,
                              number,
                              returnData,
                              savePredictions,
                              allowParallel)
    modelsRes <- list(modelsRes)
    names(modelsRes) <- classes[1]
    
    
  }else{
    modelsRes <- pblapply(classes, .trainModel,
                          object,
                          model,
                          reduction,
                          preProcess,
                          resampleMethod,
                          tuneLength,
                          seed,
                          metric,
                          number,
                          returnData,
                          savePredictions,
                          allowParallel)
    names(modelsRes) <- classes
  }
  
  cat(crayon::green("DONE!\n"))
  
  if(is.null(reclassify)){
    object@train <- modelsRes
  }else{
    object@train[names(modelsRes)] <- modelsRes
  }
  
  if(object_class == "Seurat"){
    seurat_object@misc$scPred <- object
    seurat_object
    
  }else{
    object
  }
}

.trainModel <- function(positiveClass,
                        spmodel,
                        model,
                        reduction,
                        preProcess,
                        resampleMethod,
                        tuneLength,
                        seed,
                        metric,
                        number,
                        returnData,
                        savePredictions, 
                        allowParallel){
  
  
  
  
  if(nrow(spmodel@features[[positiveClass]]) == 0){
    message("No informative principal components were identified for class: ", positiveClass)
  }
  
  names_features <- as.character(spmodel@features[[positiveClass]]$feature)
  features <- scPred:::subsetMatrix(spmodel@cell_embeddings, names_features)
  response <- as.character(spmodel@metadata$response)
  
  
  i <- response != .make_names(positiveClass)
  response[i] <- "other"
  response <- factor(response, levels = c(.make_names(positiveClass), "other"))
  
  
  if(!is.null(seed)) set.seed(seed)
  
  if(metric == "ROC"){
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           summaryFunction = twoClassSummary,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = allowParallel,
                           sampling='smote')
    
  }else if(metric == "PR"){
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           summaryFunction = prSummary,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = allowParallel,
                           sampling='smote')
    metric <- "AUC"
  }else{
    trCtrl <- trainControl(classProbs = TRUE,
                           method = resampleMethod,
                           number = number,
                           returnData = returnData,
                           savePredictions = savePredictions,
                           allowParallel = allowParallel,
                           sampling='smote')
  }
  
  
  if(metric == "AUC"){
    fit <- train(x = features,
                 y = response,
                 method = model,
                 metric = metric,
                 trControl = trCtrl,
                 preProcess = preProcess, 
                 tuneLength = tuneLength)
  }else{
    fit <- train(x = features,
                 y = response,
                 method = model,
                 preProcess= preProcess,
                 metric = metric,
                 trControl = trCtrl, 
                 tuneLength = tuneLength)
  }
  
  fit
}







#####

#mouse=qread("/rds/general/user/ems2817/home/thesis/ML/mouse_seu_filtered.qs")
mouse=qs::qread("/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/aled/mouse_seu_filtered.qs")

test=mouse[,which(mouse@meta.data$clusters=='1'|mouse@meta.data$clusters=='2'|mouse@meta.data$clusters=='3'|mouse@meta.data$clusters=='4')]


#test@meta.data$clusters=ifelse(test@meta.data$clusters=='3','DAM','non-DAM')

test<- test%>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP( dims = 1:15)

ElbowPlot(test)
DimPlot(test, group.by = "clusters", label = TRUE, repel = TRUE)

#test <- getFeatureSpace(test, "clusters")
test <-getFeatureSpace(test, "clusters", correction = "fdr", sig = 0.05, reduction = "pca")


test <- trainModel_new(test)
# mouse
test
saveRDS(test,'mouse_model_SMOTE2_aim11.rds')


get_scpred(test)
#plot_probabilities(test)

get_classifiers(test)
#caret::plot.train(get_classifiers(test)[["DAM"]])
caret::plot.train(get_classifiers(test)[["1"]])
caret::plot.train(get_classifiers(test)[["2"]])
caret::plot.train(get_classifiers(test)[["3"]])
caret::plot.train(get_classifiers(test)[["4"]])

#### cell classification

#seu <- qs::qread("/rds/general/user/ems2817/home/thesis/clara_aim12/obj_clusters_manifest_20pcs.qs")
seu <- qs::qread("/rds/general/user/ems2817/home/thesis/new_aim11/obj_clusters_20pcs.qs")
sce <- qs::qread('/rds/general/project/ukdrmultiomicsproject/live/Users/Nurun/gerrits_smith/processed_sce/Micro_sce.qs')
idx <- which(seu$barcode %in% sce$barcode)
idx <- seu$barcode[idx]
query<- subset(seu, cells = names(idx))
#
#seu=seu[,which(seu@meta.data$species=='human')]
#
DefaultAssay(query) <- "RNA"
query <- NormalizeData(query)
#
query <- scPredict(query,test)
DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")
human <- RunUMAP(query, reduction = "scpred", dims = 1:20)
#
DimPlot(human, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(query, group.by = "scpred_prediction", label = TRUE, repel = TRUE)
DimPlot(query, group.by = "integrated_snn_res.0.015", label = TRUE, repel = TRUE)
#
#
crossTab(query, "integrated_snn_res.0.015", "scpred_prediction", output = "prop")
#
qsave(query ,'/rds/general/user/ems2817/home/thesis/clara_aim12/fun/svm/query_svm2_aim11.qs')
#FeaturePlot(query, c("DAM","non-DAM"))
query@meta.data$scpred_1=signif(query@meta.data$scpred_1,2)
query@meta.data$scpred_2=signif(query@meta.data$scpred_2,2)
query@meta.data$scpred_3=signif(query@meta.data$scpred_3,2)
query@meta.data$scpred_4=signif(query@meta.data$scpred_4,2)

feature=scPred::FeaturePlot(query, c("scpred_1","scpred_2","scpred_3","scpred_4"))
ggsave("feature-svm.png", feature, height = 5, width = 7, dpi = 320)

crossTab(query, "integrated_snn_res.0.015", "scpred_prediction", output = "prop")


DimPlot(query, reduction = "umap",group.by = "scpred_1",label = TRUE,label.size = 6)
FeaturePlot(object = query, features = "scpred_1")

