## This scrip will check the top DE genes in the seasonal DGE.
## it will create plots expression level vs date for the top genes 

library(dplyr)
library(edgeR)
library(lubridate)

setwd("~/Work/RP2/ATLANTIS")

de.results.2.seasons <-  read.csv('./Season/Season_new_date/DE_genes_15Nov_winter_spring.csv')

master.Table <- read.csv("./Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv") %>%
  mutate(VISIT.DATE = as.Date(VISIT.DATE))

# expression data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE)%>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()

design <- model.matrix(~0 + new_seasons + age + gender + smoking.status + asthma.status, data = master.Table)

DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")

cpm_norm <- cpm(DGEL, normalized.lib.sizes = TRUE, log = TRUE)

## plot cpm vs date for top genes
# up genes
up_genes_cpm <- cpm_norm[c(de.results.2.seasons$Gene[(de.results.2.seasons$logFC > 0 & de.results.2.seasons$FDR < 0.05)][1:5]), ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_id") %>%
  left_join(de.results.2.seasons %>%
              dplyr::select(c(Gene, external_gene_name)), by = c("ensembl_id" = "Gene")) %>%
  tibble::column_to_rownames("external_gene_name") %>%
  dplyr::select (-ensembl_id) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample")

plot_df_up <- up_genes_cpm %>%
  tidyr::pivot_longer (c(-Sample), names_to = "Gene", values_to = "expression") %>%
  left_join(master.Table, by = c("Sample" = "GenomeScan_ID"))

plt_up <- ggplot(plot_df_up, aes(VISIT.DATE, expression)) +
  geom_point(aes(color = new_seasons, shape = asthma.status))+
  geom_smooth(colour = "black")+
  scale_colour_manual(values=c("winter_spring"= "#BFCCB5","summer_autumn"="#FFA559"),
                      aesthetics = "colour",
                      name='',
                      labels = c('summer-autumn','winter-spring')) +
  scale_shape_manual(values=c("A"= 1,"H"= 2 ),
                     aesthetics = "shape",
                     name='',
                     labels = c('Asthmatic subjects', 'Healthy subjects')) +
  ylab('Upregulated in winter-spring genes\nlog2cpm normalized expression') +
  xlab('Date') +
  facet_grid(Gene ~ ., scales="free_y")

# down genes 
down_genes_cpm <- cpm_norm[c(de.results.2.seasons$Gene[(de.results.2.seasons$logFC < 0 & de.results.2.seasons$FDR < 0.05)][1:5]), ] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ensembl_id") %>%
  left_join(de.results.2.seasons %>%
              dplyr::select(c(Gene, external_gene_name)), by = c("ensembl_id" = "Gene")) %>%
  tibble::column_to_rownames("external_gene_name") %>%
  dplyr::select (-ensembl_id) %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample")

plot_df_down <- down_genes_cpm %>%
  tidyr::pivot_longer (c(-Sample), names_to = "Gene", values_to = "expression") %>%
  left_join(master.Table, by = c("Sample" = "GenomeScan_ID"))

plt_down <- ggplot(plot_df_down, aes(VISIT.DATE, expression)) +
  geom_point(aes(color = new_seasons, shape = asthma.status))+
  geom_smooth(colour = "black")+
  scale_colour_manual(values=c("winter_spring"= "#BFCCB5","summer_autumn"="#FFA559"),
                      aesthetics = "colour",
                      name='',
                      labels = c('summer-autumn','winter-spring')) +
  scale_shape_manual(values=c("A"= 1,"H"= 2 ),
                     aesthetics = "shape",
                     name='',
                     labels = c('Asthmatic subjects', 'Healthy subjects')) +
  ylab('Downregulated in winter-spring genes\nlog2cpm normalized expression') +
  xlab('Date') +
  facet_grid(Gene ~ ., scales="free_y")

# combine up and down 
figure_all <- ggarrange(plt_up, plt_down,
                        labels = c("A", "B"),
                        common.legend = TRUE,
                        legend = "bottom",
                        ncol = 2)
png("./Season/Season_new_date/plots/top_genes_vs_date_supl3.png", width=1200, height=1500, res = 150)
print(figure_all)
dev.off()


