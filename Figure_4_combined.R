
library(ggpubr)

setwd("~/Work/RP2/ATLANTIS")

load("./Season/Season_new_date/plots/FDR_asthma_vs_asthma+season_paired_plot.Rdata")
load("./Season/Season_new_date/plots/CIBERSORT_seasons_15Nov.Rdata")

plt_narrow <- plt + 
  theme(plot.margin = unit(c(0.5, 2, 0.5, 2), "cm"))  # top, right, bottom, left
figure4_combined_gg <- ggarrange(plt_narrow,
                                 bxp_seasons, 
                                 ncol = 1,
                                 nrow = 2,
                                 heights = c(2.5, 3),
                                 labels = c("A", "B"),
                                 label.x = 1,      # Horizontal position (0 = left, 1 = right)
                                 label.y = 1,      # Vertical position (0 = bottom, 1 = top)
                                 hjust = 2.5,     # Fine-tune horizontal justification
                                 vjust = 2.5 )


png("./Season/Season_new_date/plots/Figure4_combined.png", 
    width = 1000, 
    height = 1200,
    res = 150)
print(figure4_combined_gg)
dev.off()
