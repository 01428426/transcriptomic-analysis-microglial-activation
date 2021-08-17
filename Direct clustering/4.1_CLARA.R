library(cluster)
library(factoextra)
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(argparse)
library(qs)

########
##  Parse command-line arguments                                            ####
# create parser object
parser <- ArgumentParser()
# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")
required$add_argument(
  "--cluster",
  help = "an integer",
  metavar = 1,
  required = TRUE
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()
# set working directory
cluster <- args$cluster


########
obj=qread('/rds/general/project/hda_students_data/live/Group7/General/Eleonore/thesis/new_aim11/obj_clusters_20pcs.qs')
#obj=obj[,1:500000]

# cell.embeddings: stores the coordinates for each cell in low-dimensional space.
#ifelse(dir.exists("clara_para_500samples"),"",dir.create("clara_para_500samples"))

data=data.frame(Embeddings(obj[['umap']])[,1:2])

clara.res <- clara(data, cluster, samples = 100, pamLike = TRUE, sampsize=5000)
saveRDS(clara.res,paste0("samples_umap1_500/clara_res",cluster,".rds"))

 plot=fviz_silhouette(clara.res , label = FALSE, print.summary = TRUE)


ggsave(
  plot = plot,
  filename = paste0("samples_umap1_500/clara",cluster,".png"),
  dpi = 300,
  height = 10,
  width = 10,
  units = "in")

c <- cbind(obj@meta.data$barcode, cluster = clara.res$cluster)

saveRDS(c, paste0("samples_umap1_500/cluster_para",cluster,".rds"))


