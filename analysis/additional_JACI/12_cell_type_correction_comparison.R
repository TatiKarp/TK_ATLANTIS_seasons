
library(dplyr)
library(config)

# location 
conf <- config::get()


# load data
de_seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))
de_seasons_cell_type_cor <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_season_ciliated_cor.csv"))

# df fot the plot
de_combined <- data.frame(Gene = union(de_seasons$Gene, de_seasons_cell_type_cor$Gene)) %>%
  left_join(de_seasons, by = c("Gene" = "Gene")) %>%
  left_join(de_seasons_cell_type_cor %>%
              rename_with(~ paste0(., "_corrected")), by = c("Gene" = "Gene_corrected")) %>%
  mutate(dir_FDR = ifelse(logFC > 0, 1,-1)*-log10(FDR),
         dir_FDR_corrected = ifelse(logFC_corrected > 0, 1,-1)*-log10(FDR_corrected)) %>%
  mutate(sign = case_when(
    FDR < 0.05 & FDR_corrected < 0.05 ~ "both",
    FDR < 0.05 & FDR_corrected > 0.05 ~ "initial",
    FDR > 0.05 & FDR_corrected < 0.05 ~ "corrected",
    .default = "none"))

plt <- ggplot(de_combined, aes(dir_FDR, dir_FDR_corrected)) +
  geom_point(aes(color = sign), size = 1, alpha = 0.6) +
  xlim(-5, 7.5) +
  ylim(-5, 7.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") 

png(file.path(conf$data_path, "Season/Season_new_date/plots/p_vs_p_cell_type.png"), width = 800, height = 800, res = 150)
print(plt)
dev.off()
                                     