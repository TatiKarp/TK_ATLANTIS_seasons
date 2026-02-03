library(ggpubr)
library(patchwork)
library(config)

# location 
conf <- config::get()

# load data
load(file.path(conf$data_path, "Season/Season_new_date/plots/up_PIAMA_ATLANTIS_seasons_GSVA.Rdata"))
load(file.path(conf$data_path, "Season/Season_new_date/plots/down_PIAMA_ATLANTIS_seasons_GSVA.Rdata"))
load(file.path(conf$data_path, "Season/Season_new_date/plots/Date_vs_gsva_up_down_PIAMA.Rdata"))

all <- annotate_figure(all,
                       left = text_grob("GSVA score",  rot = 90, face = "bold"))
up_down_combined <- ggarrange(
  up,
  down,         # Remove duplicate y-axis label for cleaner look
  ncol = 1,
  nrow = 2,
  common.legend = TRUE,           # Shared legend between up and down
  legend = "bottom"               # Legend position
)
up_down_combined <- annotate_figure(up_down_combined,
                                    left = text_grob("GSVA score",  rot = 90, face = "bold"))


layout_design <- "
  #A#
  ###
  BBB
"

combined_plot <- all + up_down_combined +
  plot_layout(
    design = layout_design,
    widths = c(1,4.5,1),   # Adjust to make right plots thin
    heights = c(1.6, 0.1, 1.6),  # Adjust to make top plots short
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = list(c("A", "B"))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    legend.position = "bottom"
  )


png(file.path(conf$data_path, "Season/Season_new_date/plots/Figure3_combined.png"), 
    width = 1100, 
    height = 1800,
    res = 150)
print(combined_plot)
dev.off()
