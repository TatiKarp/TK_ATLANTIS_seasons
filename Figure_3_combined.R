
library(ggpubr)

setwd("~/Work/RP2/ATLANTIS")


load("./Season/Season_new_date/plots/up_PIAMA_ATLANTIS_seasons_GSVA.Rdata")
load("./Season/Season_new_date/plots/down_PIAMA_ATLANTIS_seasons_GSVA.Rdata")
load("./Season/Season_new_date/plots/Date_vs_gsva_up_down_PIAMA.Rdata")

up_down_combined <- ggarrange(
  up,
  down,         # Remove duplicate y-axis label for cleaner look
  ncol = 1,
  nrow = 2,
  common.legend = TRUE,           # Shared legend between up and down
  legend = "right"               # Legend position
)
figure3_combined_gg <- ggarrange(all,
                                 up_down_combined, 
                                 ncol = 1,
                                 nrow = 2,
                                 heights = c(2.5, 3),
                                 labels = c("A", "B"),
                                 label.x = 1,      # Horizontal position (0 = left, 1 = right)
                                 label.y = 1,      # Vertical position (0 = bottom, 1 = top)
                                 hjust = 2.5,     # Fine-tune horizontal justification
                                 vjust = 2.5 )

annotated_figure_3 <- annotate_figure(figure3_combined_gg,
                                      left = text_grob("GSVA value",  rot = 90, face = "bold"))

png("./Season/Season_new_date/plots/Figure3_combined.png", 
    width = 1100, 
    height = 1800,
    res = 150)
print(annotated_figure_3)
dev.off()
