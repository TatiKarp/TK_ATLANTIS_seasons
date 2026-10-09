library(patchwork)
library(ggpubr)
library(config)

# location 
conf <- config::get()

# load data
annotated_figure_gsva_vs_date <- readRDS(file.path(conf$data_path, "Season/Season_new_date/plots/Date_vs_gsva_up_down_minus2samples.RDS"))

## asthma/ allergy per season
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//up_asthma_seasons_GSVA_minus2samples.Rdata"))
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//down_asthma_seasons_GSVA_minus2samples.Rdata"))
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//up_allergy_seasons_GSVA_minus2samples.Rdata"))
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//down_allergy_seasons_GSVA_minus2samples.Rdata"))


up_genes <- ggarrange(up_asthma + theme(legend.position = "none"),
                      up_allergy + rremove("ylab") + theme(legend.position = "none"),
                      common.legend = FALSE)

down_genes <- ggarrange(down_asthma, 
                        down_allergy + rremove("ylab"),
                        common.legend = TRUE,
                        legend = "bottom")

figure2_combined_gg <- ggarrange(annotated_figure_gsva_vs_date,
                                 up_genes, down_genes, 
                                 ncol = 1,
                                 nrow = 3,
                                 heights = c(2.8, 1.5, 1.5),
                                 labels = c("A", "B", ""),
                                 label.x = 1,      # Horizontal position (0 = left, 1 = right)
                                 label.y = 1,      # Vertical position (0 = bottom, 1 = top)
                                 hjust = 2.5,     # Fine-tune horizontal justification
                                 vjust = 2.5 )

annotated_figure_2 <- annotate_figure(figure2_combined_gg,
                                      left = text_grob("GSVA value",  rot = 90, face = "bold"))
png(file.path(conf$data_path,"Season/Season_new_date/plots/SF2_combined_minus2samples.png"), 
    width = 1100, 
    height = 1800,
    res = 150)
print(annotated_figure_2)
dev.off()
