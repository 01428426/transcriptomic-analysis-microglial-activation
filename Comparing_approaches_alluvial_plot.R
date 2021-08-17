library(qs)
library(Seurat)
library(SeuratObject)
library(SummarizedExperiment)
library(ggalluvial)
library(ggplot2)

seu=qread('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/obj_clusters_20pcs.qs')
seu2=qread('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/clara_aim12/obj_clusters_manifest_20pcs.qs')
seu3=qread('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/clara_aim12/query.qs')

seu@meta.data$nuclei_id <- gsub("^[a-z]{5}_", "", seu@meta.data$barcode)
seu2@meta.data$nuclei_id <- gsub("^[a-z]{5}_", "", seu2@meta.data$barcode)
seu3@meta.data$nuclei_id <- gsub("^[a-z]{5}_", "", seu3@meta.data$barcode)
cluster1= data.frame(barcode=seu@meta.data$nuclei_id, cluster1= as.numeric(seu@meta.data$integrated_snn_res.0.015))
cluster2= data.frame(barcode=seu2@meta.data$nuclei_id, cluster2= as.numeric(seu2@meta.data$integrated_snn_res.0.01))
cluster3= data.frame(barcode=seu3@meta.data$nuclei_id, cluster3= seu3@meta.data$scpred_prediction)

cluster1$cluster1=(ifelse(cluster1$cluster1=='4', 'DAM1',
                          ifelse(cluster1$cluster1=='2','DAM2', 'non-DAM')))
cluster2$cluster2=(ifelse(cluster2$cluster2=='3', 'DAM2',
                        ifelse(cluster2$cluster2=='4','DAM1', 'non-DAM')))
cluster3$cluster3=(ifelse(cluster3$cluster3=='3', 'DAM', 'non-DAM'))

alluvial=merge(cluster1, cluster2, by='barcode')
alluvial2=merge(alluvial,cluster3,by='barcode')


df=as.data.frame(table(alluvial2$cluster2,alluvial2$cluster1,alluvial2$cluster3))

ggplot(data = df,
       aes(axis1 = Var2, axis2 = Var1, axis3=Var3, y = Freq)) +
  geom_alluvium(aes(fill = Var1)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)))+
  theme_void()