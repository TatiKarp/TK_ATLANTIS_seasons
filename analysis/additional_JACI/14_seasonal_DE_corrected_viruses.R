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
 
# expression data
expression.data <- read.csv(file.path(conf$data_path, "Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header = TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID)) %>%
  as.matrix()

# viral data

# save files from CZID
files <- list.files(path = file.path(conf$data_path, "Microbes/Output/nasal_brushes_score_check"), pattern="*.csv", full.names=TRUE, recursive=FALSE)
# save samples which fits filtering conditions: 
## presence of enterovirus with z_score > 3 ##########
## presence of rhinovirus in enterovirus genus
samples <- c()
for (file in files){
  sample <- read.csv(file)
  filtered <- sample%>%
    mutate(max_z_score = pmax(new_nt_z_score, new_nr_z_score, na.rm = TRUE))%>%
    filter(max_z_score > 3)%>%
    filter(genus_tax_id == '12059') ## enterovirus 
  if (any(startsWith(filtered$name, 'Rhino'))){ #save only samples with Rhinoviruses present
    file_ID <- str_split(file, '/')[[1]][10]
    sample_ID <- str_split(file_ID, '_')[[1]][1]
    samples <- append (samples, sample_ID)
    print(filtered$name)
  }
}
## add data to master table ####
master.Table <- master.Table %>%
  mutate(Rhinovirus_presence = if_else(original_id %in% samples, 'yes', 'no'))

##### 1. DE analysis season + asthma + viral presence - main model 
design <- model.matrix(~new_seasons + age + gender + smoking.status + asthma.status + Rhinovirus_presence, data = master.Table)
edger <- run_edger_pipeline(expression.data, design)

# Conduct genewise statistical tests for a given coefficient or contrast.
# contrasts <- limma::makeContrasts(season_binary = new_seasonswinter_spring - new_seasonssummer_autumn,
#                                  levels = design)
qlf_seasons <- edgeR::glmQLFTest(edger$fit, coef = 2)
summary(decideTests(qlf_seasons))

# new_seasonswinter_spring
# Down                        722
# NotSig                    16706
# Up                          905

qlf_virus <- edgeR::glmQLFTest(edger$fit, coef = 8)
summary(decideTests(qlf_virus))

# Rhinovirus_presenceyes
# Down                        5
# NotSig                  18162
# Up                        166

# add gene names
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_ids <- rownames(edger$DGEL$counts)
all_new_gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)

de_result_seasons <- topTags(qlf_seasons, n = nrow(edger$DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id"))

de_result_viral <- topTags(qlf_virus, n = nrow(edger$DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id"))

## combine with initial season result
de_seasons_initial <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))


de_combined <- data.frame(Gene = union(de_result_seasons$Gene, de_seasons_initial$Gene)) %>%
  left_join(de_seasons_initial %>%
              rename_with(~ paste0(., "_season_init")), by = c("Gene" = "Gene_season_init")) %>%                         
  left_join(de_result_seasons %>%
              rename_with(~ paste0(., "_season")), by = c("Gene" = "Gene_season")) %>%
  left_join(de_result_viral %>%
              rename_with(~ paste0(., "_viral")), by = c("Gene" = "Gene_viral")) %>%
  mutate(dir_FDR_season_init = ifelse(logFC_season_init > 0, 1,-1)*-log10(FDR_season_init),
         dir_FDR_season = ifelse(logFC_season > 0, 1,-1)*-log10(FDR_season),
         dir_FDR_viral = ifelse(logFC_viral > 0, 1,-1)*-log10(FDR_viral)) %>%
  mutate(
    sign_season = case_when(
      FDR_season_init < 0.05 & FDR_season < 0.05 ~ "both",
      FDR_season_init < 0.05 & FDR_season > 0.05 ~ "initial",
      FDR_season_init > 0.05 & FDR_season < 0.05 ~ "corrected",
      .default = "none"),
    sign_season_virus = case_when(
      FDR_season < 0.05 & FDR_viral < 0.05 ~ "both",
      FDR_season < 0.05 & FDR_viral > 0.05 ~ "only_season",
      FDR_season > 0.05 & FDR_viral < 0.05 ~ "only_virus",
      .default = "none")) ##This none is ong- some genes were not even tested in the more complex model

# plots
# For genes that are in both viral and seasonal contrast 
plt <- ggplot(de_combined, aes(dir_FDR_viral, dir_FDR_season)) +
  geom_point(aes(color = sign_season_virus, shape = sign_season), size = 1, alpha = 0.6) +
  xlim(-5, 7.5) +
  ylim(-5, 7.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50")+
  labs(
    color = "Season & Virus sign", 
    shape = "Season sign in corrected and initial"           # renames the shape legend
  ) +
  theme(legend.position = "bottom",
        legend.box = "vertical"        # stacks legends on top of each other
  )


png(file.path(conf$data_path, "Season/Season_new_date/plots/p_vs_p_season_virus.png"), 
    width = 1000, height = 800, res = 150)
print(plt)
dev.off()

##### 2. Perform "sensitivity analysis": exclude subjects with viral reads.

# metadata
master.Table.noV <- master.Table %>%
  filter(Rhinovirus_presence == 'no')

# expression data
expression.data.noV <- expression.data[,master.Table.noV$GenomeScan_ID]

# DE analysis season
design <- model.matrix(~new_seasons + age + gender + smoking.status + asthma.status, data = master.Table.noV)
edger <- run_edger_pipeline(expression.data.noV, design)

# Conduct genewise statistical tests for a given coefficient or contrast.
# contrasts <- limma::makeContrasts(season_binary = new_seasonswinter_spring - new_seasonssummer_autumn,
#                                  levels = design)
qlf_seasons <- edgeR::glmQLFTest(edger$fit, coef = 2)
summary(decideTests(qlf_seasons))

# new_seasonswinter_spring
# Down                        296
# NotSig                    17127
# Up                          354

de_result_seasons_noV <- topTags(qlf_seasons, n = nrow(edger$DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  left_join(all_new_gene, by = c("Gene" = "ensembl_gene_id"))

## combine with initial season result
de_seasons_initial <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))


de_combined <- data.frame(Gene = union(de_result_seasons_noV$Gene, de_seasons_initial$Gene)) %>%
  left_join(de_seasons_initial %>%
              rename_with(~ paste0(., "_season_init")), by = c("Gene" = "Gene_season_init")) %>%                         
  left_join(de_result_seasons_noV %>%
              rename_with(~ paste0(., "_season_noV")), by = c("Gene" = "Gene_season_noV")) %>%
  mutate(dir_FDR_season_init = ifelse(logFC_season_init > 0, 1,-1)*-log10(FDR_season_init),
         dir_FDR_season_noV = ifelse(logFC_season_noV > 0, 1,-1)*-log10(FDR_season_noV)) %>%
  mutate(
    sign_season = case_when(
      FDR_season_init < 0.05 & FDR_season_noV < 0.05 ~ "both",
      FDR_season_init < 0.05 & FDR_season_noV > 0.05 ~ "initial",
      FDR_season_init > 0.05 & FDR_season_noV < 0.05 ~ "only_in_no_vir",
      is.na(FDR_season_init) ~ "not_tested_init",
      is.na(FDR_season_noV) ~ "not_tested_noV",
      .default = "none"))

# plots
# For genes that are in both viral and seasonal contrast 
plt <- ggplot(de_combined, aes(dir_FDR_season_init, dir_FDR_season_noV)) +
  geom_point(aes(color = sign_season), size = 1, alpha = 0.6) +
  xlim(-5, 7.5) +
  ylim(-5, 7.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50")+
  labs(
    color = "Season sign in initial and in sensitivity"          
  ) +
  theme(legend.position = "bottom")


png(file.path(conf$data_path, "Season/Season_new_date/plots/p_vs_p_season_NOvirus.png"), 
    width = 1200, height = 800, res = 150)
print(plt)
dev.off()
