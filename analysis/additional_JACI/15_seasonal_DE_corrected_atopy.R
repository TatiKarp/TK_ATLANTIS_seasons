# This code will correct the DE seasonal model for the viral presence (obtained from TNA-seq no human reads)
library(dplyr)
library(edgeR)
library(ggplot2)
library(tibble)
library(pheatmap)
library(ggrepel)
library(biomaRt)
library(config)
library(stringr)

# location 
conf <- config::get()

source(file.path(conf$script_path, "R/utils_DE.R"))

# load data
# metadata
master.Table <- read.csv(file.path(conf$data_path, "./Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv"), header = TRUE)

# atopy info
atopy_info <- read.csv(file.path(conf$data_path, "atlantis_patient_data.csv"), header =TRUE, na.strings=c("","NA")) %>%
  filter(VISIT_n == 1) %>%
  dplyr::select(c('PT','PHADRES'))

master.Table <- master.Table %>%
  left_join(atopy_info, by = c('PT'='PT')) %>%
  filter(!is.na(PHADRES))

# expression data
expression.data <- read.csv(file.path(conf$data_path, "Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header = TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID)) %>%
  as.matrix()

##### 1. DE analysis season + asthma + atopy
design <- model.matrix(~new_seasons + age + gender + smoking.status + asthma.status + PHADRES, data = master.Table)
edger <- run_edger_pipeline(expression.data, design)

# Conduct genewise statistical tests for a given coefficient or contrast.
# contrasts <- limma::makeContrasts(season_binary = new_seasonswinter_spring - new_seasonssummer_autumn,
#                                  levels = design)
qlf_seasons <- edgeR::glmQLFTest(edger$fit, coef = 2)
summary(decideTests(qlf_seasons))

# new_seasonswinter_spring
# Down                        614
# NotSig                    16913
# Up                          517

# add gene names
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_ids <- rownames(edger$DGEL$counts)
all_new_gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)

de_result_seasons_atopy_corrected <- topTags(qlf_seasons, n = nrow(edger$DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id"))

## 2. lot comparison 
de_seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))

#
de_combined <- data.frame(Gene = union(de_seasons$Gene, de_result_seasons_atopy_corrected$Gene)) %>%
  left_join(de_seasons, by = c("Gene" = "Gene")) %>%
  left_join(de_result_seasons_atopy_corrected %>%
              rename_with(~ paste0(., "_corrected")), by = c("Gene" = "Gene_corrected")) %>%
  mutate(dir_FDR = ifelse(logFC > 0, 1,-1)*-log10(FDR),
         dir_FDR_corrected = ifelse(logFC_corrected > 0, 1,-1)*-log10(FDR_corrected)) %>%
  mutate( 
    tested = case_when(
      is.na(FDR) ~ "not tested initial", 
      is.na(FDR_corrected) ~ "not tested corrected", 
      !(is.na(FDR)) & !is.na(FDR_corrected) ~ "tested both"),
    significant = case_when(
    FDR < 0.05 & FDR_corrected < 0.05 ~ "both",
    FDR < 0.05 & FDR_corrected > 0.05 ~ "initial",
    FDR > 0.05 & FDR_corrected < 0.05 ~ "corrected",
    is.na(FDR) ~ "not tested initial", 
    is.na(FDR_corrected) ~ "not tested corrected", 
    .default = "none"))

table(de_combined$significant)

plt <- ggplot(de_combined, aes(dir_FDR, dir_FDR_corrected)) +
  geom_point(aes(color = significant), size = 1, alpha = 0.6) +
  xlim(-5, 7.5) +
  ylim(-5, 7.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") 

png(file.path(conf$data_path, "Season/Season_new_date/plots/p_vs_p_atopy_corrected.png"), width = 800, height = 800, res = 150)
print(plt)
dev.off()

