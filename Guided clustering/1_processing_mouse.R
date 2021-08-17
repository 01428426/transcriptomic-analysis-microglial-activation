
### MOUSE DATASET

#Convert the mouse gene symbol to its human ortholog using the lookup table.
library(qs)
library(SingleCellExperiment)
sce_mouse=qread('sce_keren_shaul.qs')
SummarizedExperiment::rowData(sce_mouse)$human_gene_ortholog 

orthologs$Mouse_gene_name=orthologs$`Mouse gene name`
orthologs2=orthologs[!duplicated(orthologs$Mouse_gene_name), ]


vector = data.frame(Mouse_gene_name= rowData(sce_mouse)$gene)
vector$Mouse_gene_name=as.character(vector$Mouse_gene_name)
merge= left_join(vector,orthologs2, by='Mouse_gene_name')
merge=merge$`Gene name`

gene=merge$`Gene name`

rownames(sce_mouse) <- gene

#remove NA --> no human equivalent 
sce_mouse=sce_mouse[!is.na(rownames(sce_mouse)),]

rowData(sce_mouse)$gene= rownames(sce_mouse)

sce_mouse<- sce_mouse[!duplicated(SummarizedExperiment::rowData(sce_mouse)$gene), ]


sce_mouse$species <- 'mouse'

#get rid of cluster 0
sce_mouse=sce_mouse[,which(sce_mouse$clusters!=0)]
sce_mouse=sce_mouse[,which('1'|sce_mouse$cluster=='2'|sce_mouse$cluster=='3'|sce_mouse$cluster=='4')]

seu_mouse<- CreateSeuratObject(counts = SingleCellExperiment::counts(sce_mouse), 
                               meta.data = as.data.frame(SingleCellExperiment::colData(sce_mouse), project = "mouse"))

seu_mouse[["RNA"]] <- AddMetaData(seu_mouse[["RNA"]], as.data.frame(SummarizedExperiment::rowData(sce_mouse)),col.name = NULL)

qs::qsave(seu_mouse,"/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/round_2/mouse_seu.qs")
