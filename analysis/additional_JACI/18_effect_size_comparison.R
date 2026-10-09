
library(duckplyr)
library(grid)
library(ggplot2)
library(config)

# location 
conf <- config::get()

######## Compare p-values ###########

# Create one table with DE statistics out of 2 
# Arguments: 2 tables resulting from DE analysis

create_compare_df <- function(df_one, df_two, name_one, name_two) {
  overlap_genes <- intersect(df_one$Gene, df_two$Gene)
  comparison_df <- df_one %>%
    filter(Gene %in% overlap_genes) %>%
    dplyr::select(c(Gene, logFC, PValue, FDR)) %>%
    rename_with(~ paste0(name_one, "_", .x), all_of(c("logFC", "PValue", "FDR"))) %>%
    left_join(df_two %>%
                dplyr::select(c(Gene, logFC, PValue, FDR)), by = c("Gene" = "Gene")) %>%
    rename_with(~ paste0(name_two, "_", .x), all_of(c("logFC", "PValue", "FDR"))) %>%
    mutate(
      significant = case_when(
        .data[[paste0(name_one, "_FDR")]] < 0.05 & .data[[paste0(name_two, "_FDR")]] < 0.05 ~ "both models",
        .data[[paste0(name_one, "_FDR")]] < 0.05 & .data[[paste0(name_two, "_FDR")]] >= 0.05 ~ paste(name_one, "in model"),
        .data[[paste0(name_one, "_FDR")]] >= 0.05 & .data[[paste0(name_two, "_FDR")]] < 0.05 ~ paste(name_two, "in model"),
        TRUE ~ "none"
      )
    )
  
  return(comparison_df)
}

## 1. Compare models with and without seasons 

de.result.asthma.seasons <- read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_asthma.csv"))
de.result.asthma <- read.csv(file.path(conf$data_path, "Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv"))

comparison_df <- create_compare_df(de.result.asthma, de.result.asthma.seasons, "asthma", "season_asthma")


comparison_df <- comparison_df %>%
  mutate(legend_names = case_when(
    significant == "asthma in model" ~ "not accounting for seasons",
    significant == "both models" ~ "both models",
    significant == "season_asthma in model" ~ "accounting for seasons",
    significant == "none" ~ "none"
  ))

plt <- ggplot(data = comparison_df %>%
                filter(season_asthma_FDR < 0.05 | asthma_FDR<0.05), 
              aes(x=asthma_logFC, y=season_asthma_logFC)) +
  geom_point(aes(color = legend_names), size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("#B4251A", "darkgray", "#163A7D")) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", linewidth = 0.1) +
  xlab(expression("logFC, DE genes in asthma, not accounting for season")) +
  ylab(expression("logFC, DE genes in asthma, accounting for season")) +
  labs(color = "Significant DE genes in the model:") +
  theme_bw() +
  theme(legend.position="right", 
        #legend.title=element_blank(), 
        legend.text = element_text(size = rel(1)), 
        axis.title.x = element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)))

png(file.path(conf$data_path, "Season/Season_new_date/plots/logfc_asthma_vs_asthma+season.png"), 
    width = 1200, height = 800, res = 150)
print(plt)
dev.off()

# check the logFC difference and correlation 
comparison_df <- comparison_df %>%
  mutate(abs_logfc_difference = abs(season_asthma_logFC - asthma_logFC))

comparison_df_only_sign <- comparison_df %>%
  filter(legend_names != "none") 

comparison_df_only_sign%>%
  pull(abs_logfc_difference) %>%
  hist()
