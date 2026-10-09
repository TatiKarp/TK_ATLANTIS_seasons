## This script will check DEG between seasons (15th Nov) when corrected for hospital 
library(dplyr)
library(edgeR)
library(ggplot2)
library(tibble)
library(pheatmap)
library(ggrepel)

# location 
conf <- config::get()

source(file.path(conf$script_path, "R/utils_DE.R"))

# load data
# metadata
master.Table <- read.csv(file.path(conf$data_path, "./Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv"), header = TRUE)

# expression data
expression.data <- read.csv(file.path(conf$data_path, "Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header = TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID)) %>%
  as.matrix()

# hospital info
master.Table.hospital <- read.csv(file.path(conf$data_path, "patient_data_final.csv"), header =TRUE, na.strings=c("","NA")) %>%
  as.data.frame()%>%
  dplyr::mutate(
    GenomeScan_ID = stringr::str_trim(GenomeScan_ID), #remove whitespace from start/end of the string
    GenomeScan_ID = paste0("X", GenomeScan_ID),
    GenomeScan_ID = gsub("-", ".", GenomeScan_ID))

master.Table <- master.Table%>%
  left_join(master.Table.hospital%>%
              dplyr::select(c('GenomeScan_ID', 'hospital')))%>%
  mutate(hospital = as.factor(hospital))

## differential expression seasons 
design <- model.matrix(~0 + new_seasons + age + gender + smoking.status + asthma.status + hospital, data = master.Table)
edger <- run_edger_pipeline(expression.data, design)

# DGEL <- edgeR::DGEList(expression.data)
# keep <- edgeR::filterByExpr(DGEL, design) 
# DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
# DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
# ATLANTIS_logcpm_horm_hospital<- cpm(DGEL, normalized.lib.sizes=TRUE, log=TRUE)
# 
# DGEL <- edgeR::estimateDisp(DGEL, design, legacy = TRUE)
# fit <- edgeR::glmQLFit(DGEL, design)
# 
contrasts <- limma::makeContrasts(season_binary = new_seasonswinter_spring - new_seasonssummer_autumn,
                                  levels = design)

qlf <- edgeR::glmQLFTest(edger$fit, contrast = contrasts[,"season_binary"])
summary(decideTests(qlf))

# -1*new_seasonssummer_autumn 1*new_seasonswinter_spring
# Down                                                      693
# NotSig                                                  22629
# Up                                                        517

# add gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_ids <- rownames(edger$DGEL$counts)
all_new_gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)

de.results.corrected.hospital <- topTags(qlf, n = nrow(edger$DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id"))

write.csv(de.results.corrected.hospital, file.path(conf$data_path, "Season/Season_new_date/DE_genes_new_seasons_corrected_hospital.csv"), row.names = FALSE)

## compare genes corrected and not corrected for hoslital ##
de_result_seasons_hospital_corrected  <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_new_seasons_corrected_hospital.csv"))
de_seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))

de_combined <-  data.frame(Gene = union(de_seasons$Gene, de_result_seasons_hospital_corrected$Gene)) %>%
  left_join(de_seasons, by = c("Gene" = "Gene")) %>%
  left_join(de_result_seasons_hospital_corrected %>%
              rename_with(~ paste0(., "_corrected")), by = c("Gene" = "Gene_corrected")) %>%
  mutate(
    dir_FDR = case_when(
      is.na(logFC)  ~ 0,
      logFC < 0   ~ -log10(FDR),
      logFC > 0   ~ log10(FDR)),
    dir_FDR_corrected = case_when(
      is.na(logFC_corrected)  ~ 0,
      logFC_corrected < 0   ~ -log10(FDR_corrected),
      logFC_corrected > 0   ~ log10(FDR_corrected))) %>%
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

png(file.path(conf$data_path, "Season/Season_new_date/plots/p_vs_p_hospital_corrected.png"), width = 800, height = 800, res = 150)
print(plt)
dev.off()

### check for coloblindness
library(colorBlindness)
cvdPlot(hosp_correcred_comp)