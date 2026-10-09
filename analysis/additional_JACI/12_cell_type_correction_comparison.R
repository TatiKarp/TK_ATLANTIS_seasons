# this script will perfom DGE analysis - contratin. season, correcting the the ciliated cells
library(dplyr)
library(config)
library(ggplot2)

# location 
conf <- config::get()


# load data
de_seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))
de_seasons_cell_type_cor <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_season_ciliated_basal_cor.csv"))
de_cilated_cells <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_season_ciliated_basal_cor_cilContr.csv"))
de_basal_cells <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma_season_ciliated_basal_cor_basalContr.csv"))

# df fot the plot
de_combined <- data.frame(Gene = union(de_seasons$Gene, de_seasons_cell_type_cor$Gene)) %>%
  left_join(de_seasons, by = c("Gene" = "Gene")) %>%
  left_join(de_seasons_cell_type_cor %>%
              rename_with(~ paste0(., "_corrected")), by = c("Gene" = "Gene_corrected")) %>%
  left_join(de_cilated_cells %>%
              rename_with(~ paste0(., "_ciliated_cor")), by = c("Gene" = "Gene_ciliated_cor")) %>%
  left_join(de_basal_cells %>%
              rename_with(~ paste0(., "_basal_cor")), by = c("Gene" = "Gene_basal_cor")) %>%
  mutate(dir_FDR = ifelse(logFC > 0, 1,-1)*-log10(FDR),
         dir_FDR_corrected = ifelse(logFC_corrected > 0, 1,-1)*-log10(FDR_corrected)) %>%
  mutate(significant = case_when(
    FDR < 0.05 & FDR_corrected < 0.05 ~ "both",
    FDR < 0.05 & FDR_corrected > 0.05 ~ "initial",
    FDR > 0.05 & FDR_corrected < 0.05 ~ "corrected",
    is.na(FDR) ~ "not tested initial", 
    is.na(FDR_corrected) ~ "not tested corrected", 
    .default = "none"))

table(de_combined$significant)

# fdr vs fdr
plt <- ggplot(de_combined, aes(dir_FDR, dir_FDR_corrected)) +
  geom_point(aes(color = significant), size = 1, alpha = 0.6) +
  xlim(-5, 7.5) +
  ylim(-5, 7.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") 

png(file.path(conf$data_path, "Season/Season_new_date/plots/p_vs_p_cell_type.png"), width = 800, height = 800, res = 150)
print(plt)
dev.off()

# logFC vs logFC
plt <- ggplot(de_combined%>%
                filter(FDR < 0.05), aes(logFC, logFC_corrected)) +
  geom_point(aes(color = sign), size = 1, alpha = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") 

png(file.path(conf$data_path, "Season/Season_new_date/plots/logfc_vs_logfc_cell_type.png"), width = 800, height = 800, res = 150)
print(plt)
dev.off()

cor(de_combined[de_combined$FDR < 0.05, "logFC"], de_combined[de_combined$FDR < 0.05,"logFC_corrected"])
# 0.9517933
                                     