library(patchwork)
library(ggpubr)
library(config)

# location 
conf <- config::get()

# load data
annotated_figure_gsva_vs_date <- readRDS(file.path(conf$data_path, "Season/Season_new_date/plots/Date_vs_gsva_up_down.RDS"))

## asthma/ allergy per season
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//up_asthma_seasons_GSVA.Rdata"))
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//down_asthma_seasons_GSVA.Rdata"))
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//up_allergy_seasons_GSVA.Rdata"))
load(file = file.path(conf$data_path, "Season/Season_new_date/plots//down_allergy_seasons_GSVA.Rdata"))

up_genes <- ggarrange(up_asthma + theme(legend.position = "none"),
                      up_allergy + rremove("ylab") + theme(legend.position = "none"),
                      common.legend = FALSE)

down_genes <- ggarrange(down_asthma, 
                        down_allergy + rremove("ylab"),
                        common.legend = TRUE,
                        legend = "bottom")

layout_design <- "
  A
  B
" 
up_down_stack <- up_genes + down_genes +
c


gsva_label <- ggplot() +
  annotate(
    "text",
    x = 0.5, y = 0.5,
    label = "GSVA score",
    angle = 90,
    fontface = "bold",
    size = 4
  ) +
  theme_void()

up_down_wrapped <- gsva_label + up_down_stack +
  plot_layout(
    widths = c(0.08, 1)  # adjust spacing if needed
  )

## combine with patchwork:

layout_design <- "
  #A#
  ###
  BBB
"

combined_plot <- annotated_figure_gsva_vs_date + up_down_wrapped +
  plot_layout(
    design = layout_design,
    widths = c(1,4.5,1),   # Adjust to make right plots thin
    heights = c(1.6, 0.1, 1.8),  # Adjust to make top plots short
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = list(c("A", "B"))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    legend.position = "bottom"
  )

png(file.path(conf$data_path, "Season/Season_new_date/plots/Figure2_combined.png"), 
    width = 1100, 
    height = 2000,
    res = 150)
print(combined_plot)
dev.off()
