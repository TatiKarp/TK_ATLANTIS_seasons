# This script will compare identified asthma genes (before and after corection for seasons) with the genes presentedd in the metaanalysis paper:
# Tsai YH et al, Meta-analysis of airway epithelium gene expression in asthma. Eur Respir J. 2018 May 17;51(5):1701962. 
# doi: 10.1183/13993003.01962-2017.
# From this paper ^ the supplementary table 3 was used. 
# and meta analys significance was defined as  q-values < 0.05  (as in the original paper from Tsai YH et al)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(openxlsx)
library(fgsea)

# location 
conf <- config::get()

# 1.1 load reasults from the metaanalysis - only significant genes
meta_genes <- read.xlsx(file.path(conf$data_path, "/Season/Season_new_date/tsai_erj_paper/FullMetaAnalysisOutput.xlsx"),
                        colNames = T) %>%
  filter(q.value <0.05)

# load de asthma genes with and without season correction
de.result.asthma <- read.csv(file.path(conf$data_path, "Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv")) %>%
  filter(FDR < 0.05)
de.result.asthma.seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma.csv")) %>%
  filter(FDR < 0.05)

# compare asthma genes (no season correction)
overlap.a.up <- length(intersect(meta_genes[meta_genes$Hedges.g >0, "Gene.Symbol"], de.result.asthma[de.result.asthma$logFC >0, "hgnc_symbol"]))
overlap.a.s.up <- length(intersect(meta_genes[meta_genes$Hedges.g >0, "Gene.Symbol"], de.result.asthma.seasons[, "external_gene_name"]))

overlap.a.down <- length(intersect(meta_genes[meta_genes$Hedges.g <0, "Gene.Symbol"], de.result.asthma[de.result.asthma$logFC >0, "hgnc_symbol"]))
overlap.a.s.down <- length(intersect(meta_genes[meta_genes$Hedges.g <0, "Gene.Symbol"], de.result.asthma.seasons[, "external_gene_name"]))

# handle multiple gene names properly:
meta_genes <- meta_genes %>%
  mutate(group_id = row_number()) # record to track back to original gene names 

meta_genes_names <- meta_genes %>%
  separate_longer_delim("Gene.Symbol", delim = "///") %>%
  mutate(gene_match_no_s = `Gene.Symbol` %in% de.result.asthma$hgnc_symbol,
         gene_match_with_s = `Gene.Symbol` %in% de.result.asthma.seasons$external_gene_name)

result <- meta_genes_names %>%
  group_by(group_id) %>%
  summarise(
    overlap_no_s = any(gene_match_no_s),                          # TRUE/FALSE
    matched_genes_no_s = paste(`Gene.Symbol`[gene_match_no_s], collapse = "; "),  # which genes matched
    overlap_with_s = any(gene_match_with_s),
    matched_genes_with_s = paste(`Gene.Symbol`[gene_match_with_s], collapse = "; "),  # which genes matched
    .groups = "drop"
  )

#  Join back to original df, and compare the direction 
meta_genes_compared <- meta_genes %>%
  left_join(result, by = "group_id") %>%
  mutate(matched_genes_no_s = if_else(matched_genes_no_s == "", NA, matched_genes_no_s),
        matched_genes_with_s = if_else(matched_genes_with_s == "", NA, matched_genes_with_s)) %>%
  left_join(de.result.asthma, by = c("matched_genes_no_s"="hgnc_symbol"), na_matches = "never") %>%
  left_join(de.result.asthma.seasons, by = c("matched_genes_with_s" = "external_gene_name"), na_matches = "never") %>%
  mutate(direction_no_s = if_else(logFC.x * Hedges.g >0, "same", "different"),
         direction_with_s = if_else(logFC.y * Hedges.g >0, "same", "different")) # check if the direction of difference is the same 


# test: 
## chi-square 
table(master.Table.ATLANTIS$new_seasons, master.Table.ATLANTIS$Rhinovirus_presence)
chisq.test(master.Table.ATLANTIS$new_seasons, master.Table.ATLANTIS$Rhinovirus_presence, correct=FALSE)

# 1.2 add GSEA
full_list <- read.xlsx(file.path(conf$data_path, "/Season/Season_new_date/tsai_erj_paper/FullMetaAnalysisOutput.xlsx"),
                       colNames = T) %>%
  mutate(rank_column = if_else(Hedges.g > 0, -log10(q.value), log10(q.value))) %>%
  arrange(desc(rank_column))

ranks <- full_list$rank_column
names(ranks) <- full_list$Gene.Symbol
de_results_only_in_seasons <- de.result.asthma.seasons[!(de.result.asthma.seasons$Gene  %in% de.result.asthma$Gene),] %>%
  filter(!external_gene_name == "")
gene_sets <- list(Up = de_results_only_in_seasons[de_results_only_in_seasons$logFC > 0, "external_gene_name"],
                  Down = de_results_only_in_seasons[de_results_only_in_seasons$logFC < 0, "external_gene_name"])

# two tailed GSEA
fgseaRes_both <-fgsea(pathways = gene_sets,
                      stats = ranks,
                      minSize=15,
                      maxSize=2000,
                      nPermSimple = 10000)

# list  of genes from a different paper: 
#Pandey, G., Pandey, O.P., Rogers, A.J. et al. A Nasal Brush-based Classifier of Asthma Identified by Machine Learning Analysis of Nasal RNA Sequence Data. Sci Rep 8, 8826 (2018).
#https://doi.org/10.1038/s41598-018-27189-4
de_Pandey_up <- read.xlsx(file.path(conf$data_path,"/Season/Season_new_date/tsai_erj_paper/41598_2018_27189_MOESM2_ESM.xlsx"),
                          sheet = 1,
                          startRow = 4,
                          colNames = T) 
overlap.a.up <- length(intersect(de_Pandey_up[de_Pandey_up$FDR < 0.05, "Gene"], de.result.asthma[de.result.asthma$logFC >0, "hgnc_symbol"]))

overlap.a.s.up <- length(intersect(de_Pandey_up[de_Pandey_up$FDR < 0.05, "Gene"], de.result.asthma.seasons[de.result.asthma$logFC >0, "external_gene_name"]))

de_Pandey_down <- read.xlsx(file.path(conf$data_path,"/Season/Season_new_date/tsai_erj_paper/41598_2018_27189_MOESM2_ESM.xlsx"),
                            sheet = 3,
                            startRow = 4,
                            colNames = T) 

overlap.a.down <- length(intersect(de_Pandey_up[de_Pandey_up$FDR < 0.05, "Gene"], de.result.asthma[de.result.asthma$logFC <0, "hgnc_symbol"]))
overlap.a.s.down <- length(intersect(de_Pandey_up[de_Pandey_up$FDR < 0.05, "Gene"], de.result.asthma.seasons[de.result.asthma$logFC <0, "external_gene_name"]))

mat <- matrix(c(overlap.a.up+overlap.a.down, 126,
                overlap.a.s.up + overlap.a.s.down, 271), nrow = 2,
              dimnames = list(c("Overlap", "ATLANTIS"),
                              c("no_season", "with_season")))
chisq.test(mat)

# define geneSets: 
geneSets_high <- list(Up_genes = de_results_only_in_seasons[de_results_only_in_seasons$logFC>0,]$external_gene_name,
                      Down_genes = de_results_only_in_seasons[de_results_only_in_seasons$logFC<0,]$external_gene_name)

# define background
meta_genes_all <- read.xlsx(file.path(conf$data_path, "/Season/Season_new_date/tsai_erj_paper/NIHMS1607987-supplement-Supplementary_Table_3.xlsx"),
                            startRow = 2,
                            colNames = T)

ordered.de_Pandey_up <- de_Pandey_up %>%
  mutate(`-10logFDR` = -log10(FDR)) %>%
  dplyr::arrange(desc(`-10logFDR`)) 
