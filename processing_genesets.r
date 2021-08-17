library(dplyr)
setwd("~/nfancy/projects/gerrits_lau_smith")

human_to_mouse <- read.delim("~/nfancy/projects/ad_unsorted/52_samples/genesets/mart_expor_human_mouse_orthologt.txt")

###
keren_shaul <- readxl::read_xlsx("genesets/Keren-shaul_2017/mmc3.xlsx")
keren_shaul <- keren_shaul %>%
  mutate(`Fold-change (DAM to homeostatic microglia)` = as.numeric(`Fold-change (DAM to homeostatic microglia)`),
         `-log10(DAM) p-value (Mann-Whitney)` = as.numeric(`-log10(DAM) p-value (Mann-Whitney)`),
         `DAM FDR p-value` = as.numeric(`DAM FDR p-value`)) %>%
  filter(!is.nan(`DAM FDR p-value`) )%>%
  filter(!is.nan(`Fold-change (DAM to homeostatic microglia)`)) %>%
  mutate(logFC = as.numeric(`Fold-change (DAM to homeostatic microglia)`),
         pvalue= 10^(-1*`-log10(DAM) p-value (Mann-Whitney)`),
         padj = p.adjust(pvalue)) %>%
  filter(padj <= 0.1) %>%
  #filter(logFC > 0) %>%
  dplyr::rename(Mouse.gene.name = `Gene name`)

keren_shaul <- left_join(keren_shaul, human_to_mouse)
DAM <- keren_shaul %>%
  filter(!is.na(Gene.name)) %>%
  dplyr::rename(gene_name = Gene.name) %>%
  select(c(gene_name, logFC, pvalue, padj)) %>%
  filter(logFC > 0)

Homeostatic <- keren_shaul %>%
  filter(!is.na(Gene.name)) %>%
  dplyr::rename(gene_name = Gene.name) %>%
  select(c(gene_name, logFC, pvalue, padj)) %>%
  filter(logFC < 0)

# keren_shaul_1 <- readxl::read_xlsx("genesets/Keren-shaul_2017/mmc2.xlsx")
# Micro1 <- keren_shaul_1 %>%
#   mutate(`-log10(p-value) (Microglia1 vs. Microglia3)` = as.numeric(`-log10(p-value) (Microglia1 vs. Microglia3)`)) %>%
#   filter(!is.nan(`-log10(p-value) (Microglia1 vs. Microglia3)`)) %>%
#   filter(`up/down` == 1) %>%
#   mutate(pvalue= 10^(-1*`-log10(p-value) (Microglia1 vs. Microglia3)`),
#          padj = p.adjust(pvalue)) %>%
#   dplyr::rename(Mouse.gene.name = `...1`)
#   
# Micro3 <- keren_shaul_1 %>%
#   mutate(`-log10(p-value) (Microglia1 vs. Microglia3)` = as.numeric(`-log10(p-value) (Microglia1 vs. Microglia3)`)) %>%
#   filter(!is.nan(`-log10(p-value) (Microglia1 vs. Microglia3)`)) %>%
#   filter(`up/down` == -1) %>%
#   mutate(pvalue= 10^(-1*`-log10(p-value) (Microglia1 vs. Microglia3)`),
#          padj = p.adjust(pvalue)) %>%
#   dplyr::rename(Mouse.gene.name = `...1`)
  

##
Olah <- readxl::read_xlsx("genesets/Olah_2018/41467_2018_2926_MOESM6_ESM.xlsx")
Aging <- Olah %>%
  filter(adj.P.Val <= 0.001, logFC >= 1) %>%
  dplyr::rename(gene_name = ID,
                pvalue = P.Value,
                padj = adj.P.Val) %>%
  select(c(gene_name, logFC, pvalue, padj)) 


##
Galatro_table_1 <- readxl::read_xlsx("genesets/Galatro_2017/41593_2017_BFnn4597_MOESM4_ESM.xlsx",
                                     skip = 1)

Galatro_table_1 <- Galatro_table_1 %>%
  select(c(external_gene_name, Median_Brain, 
           Median_microglia, `glia-brain_logFC`, `P.Value...73`,
           `adj.P.Val...74`))

Core_Micro1 <- Galatro_table_1 %>%
  select(c(external_gene_name, Median_Brain, 
  Median_microglia, `glia-brain_logFC`, `P.Value...73`,
  `adj.P.Val...74`)) %>%
  dplyr::rename(gene_name = external_gene_name,
                logFC = `glia-brain_logFC`,
                pvalue =  P.Value...73,
                padj = adj.P.Val...74) %>%
  filter(padj <= 0.001, logFC <= 4) %>%
  select(c(gene_name, logFC, pvalue, padj)) %>%
  arrange(padj)


Galatro_table_2 <- readxl::read_xlsx("genesets/Galatro_2017/41593_2017_BFnn4597_MOESM7_ESM.xlsx",
                             sheet = "microglia-macrophage", skip = 1)
Macrophage <- Galatro_table_2 %>%
  filter(adj.P.Val <= 0.001, logFC <= -2) %>%
  dplyr::rename(gene_name = external_gene_name,
                pvalue = P.Value,
                padj = adj.P.Val) %>%
  select(c(gene_name, logFC, pvalue, padj)) %>%
  arrange(padj)


Core_Micro2 <- Galatro_table_2 %>%
  filter(adj.P.Val <= 0.001, logFC >= 2) %>%
  dplyr::rename(gene_name = external_gene_name,
                pvalue = P.Value,
                padj = adj.P.Val) %>%
  select(c(gene_name, logFC, pvalue, padj)) %>%
  arrange(padj)

##

Patir <- readxl::read_xlsx("genesets/Patir_2018/glia23572-sup-0001-tables3.xlsx",
                           skip = 3)
Core_Micro3 <- Patir %>%
  filter(`Core signature/ overlap in two dataset` == "Core") %>%
  dplyr::rename(gene_name = "Gene symbol") %>%
  select(gene_name)

save.image("genesets/curated_genesets.rda")

