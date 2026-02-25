
library(ggpubr)
library(patchwork)
library(config)

# location 
conf <- config::get()

# load data
load(file.path(conf$data_path, "Season/Season_new_date/plots/FDR_asthma_vs_asthma+season_paired_plot.Rdata"))
load(file.path(conf$data_path, "Season/Season_new_date/plots/CIBERSORT_seasons_15Nov.Rdata"))

plt_nolegend <- plt + guides(fill = "none", color = "none", linetype = "none")


layout_design <- "
  #A#
  BBB
" 

combined_plot <- plt_nolegend + bxp_seasons +
  plot_layout(
    design = layout_design,
    heights = c(2.25, 2.5),
    widths = c(0.1,2,0.1) # Adjust to make top plots short
    ) +
  plot_annotation(
    tag_levels = list(c("A", "B"))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 16),
    legend.position = "bottom"
  )


png(file.path(conf$data_path, "Season/Season_new_date/plots/Figure4_combined.png"), 
    width = 1000, 
    height = 1200,
    res = 150)
print(combined_plot)
dev.off()
