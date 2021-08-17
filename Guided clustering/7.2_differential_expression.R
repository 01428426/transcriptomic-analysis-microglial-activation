#!/usr/bin/env Rscript
# Perform differential expression analysis 
# Nurun Fancy <n.fancy@imperial.ac.uk>
#   ____________________________________________________________________________
#   Initialization                                                          ####
##  ............................................................................
##  Load packages                                                           ####
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(purrr)
library(argparse)
library(scFlow)
library(qs)
##  ............................................................................
##  Parse command-line arguments                                            ####
# create parser object
parser <- ArgumentParser()
# specify options
required <- parser$add_argument_group("Required", "required arguments")
optional <- parser$add_argument_group("Optional", "required arguments")
required$add_argument(
  "--sce",
  help = "path to sce.qs object",
  metavar = "/dir/sce/sce.qs",
  required = TRUE
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
args <- parser$parse_args()
# set working directory
setwd("/data/thesis/clara_aim12")
ensembl_fp <- "/data/thesis/ensembl_mappings_human.tsv" 
dir.create("de_new_aim12")
sce <- qs::qread(args$sce)
#starting differential expression on variable
contrast <- "diagnosis"
#mod1 <- group+(1|individual)+cngeneson+pc_mito + brain_region(when all brains are together) + sex + RIN
#mod2 <- aBeta+(1|individual)+cngeneson+pc_mito
#mod3 <- aBeta+(1|individual)+cngeneson+pc_mito+pTau+sex
mod <- "glmer_mod1"
celltype <- gsub("_sce", "", tools::file_path_sans_ext(basename(args$sce)))
outdir_name <- sprintf("%s/%s_%s/%s", "de", contrast, mod, celltype)
dir.create(outdir_name, recursive = TRUE)
confounders <- c("cngeneson",
                 "pc_mito",
                 "brain_region",
                 "sex",
                 "dataset")
de_results <- perform_de(
  sce = sce,
  de_method = "MASTZLM",
  ebayes = FALSE,
  mast_method = "glmer",
  min_counts = 1,
  min_cells_pc = 0.10,
  rescale_numerics = TRUE,
  dependent_var = contrast,
  ref_class = "CTR",
  confounding_vars = confounders,
  random_effects_var = "manifest",
  ensembl_mapping_file = ensembl_fp,
  unique_id_var = "individual",
  subset_var = "diagnosis",
  subset_class = "CTR",
  nAGQ = 0
)
for (result in names(de_results)) {
  file_name <- sprintf("%s.tsv", result)
  write.table(de_results[[result]],
              file = file.path(outdir_name, file_name),
              sep = "\t", quote = FALSE, 
              col.names = TRUE, row.names = FALSE)
}



