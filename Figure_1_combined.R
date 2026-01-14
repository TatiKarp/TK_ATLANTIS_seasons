
# This script will create figure one for paper:
# combinig selection of the season and volcanoplot 

library(patchwork)
library(lubridate)

setwd("~/Work/RP2/ATLANTIS")

master.Table <- read.csv("./Season/Season_new_date/ATLANTIS_master_table_seasons_new.csv")

## define first day of each month 
st <- as.Date("2015-01-01")
en <- as.Date("2015-12-01")
first_month <- seq(st + 1, en + 1, by = "1 month") - 1
yday(first_month)

## add 15th of each month
breaks_df <- c()
for (day_n in first_month) {
  print(yday(as.Date(day_n)))
  breaks_df <- c(breaks_df, yday(as.Date(day_n)), yday(as.Date(day_n)) + 14)
}

# Load DE result with 15 days shift (from code ~/ATLANTIS_project/Season/Season_new_date/seasons_DE_15_shift.R)
load("./Season/Season_new_date/15days_shift.Rdata")

dates_genes_df <- data.frame(Date = as.Date(breaks_df[1:length(n_deg)]-1, origin = "2015-01-01"),
                             n_gene = results_df$n_total, 
                             n_up = results_df$n_up,
                             n_down = results_df$n_down) %>%
  mutate(Date_plus6months = Date %m+% months(6)) %>%
  mutate(Date_start = paste0(format(Date, "%d%b"))) %>%
  mutate(Date_start = factor(Date_start, levels = Date_start)) %>%
  mutate(Date_range = paste0(format(Date, "%d%b"), "-", format(Date_plus6months, "%d%b"))) %>%
  mutate(Date_range = factor(Date_range, levels = Date_range))

# Plot
figure1 <- ggplot(dates_genes_df, aes(x = Date_start)) +
  geom_line(aes(y = n_up, color = "Higher expressed"), size = 0.8, group = 1) +
  geom_point(aes(y = n_up, color = "Higher expressed"), size = 2) +
  
  geom_line(aes(y = n_down, color = "Lower expressed"), size = 0.8, group = 1) +
  geom_point(aes(y = n_down, color = "Lower expressed"), size = 2) +
  
  geom_line(aes(y = n_gene, color = "Total DEGs"), size = 0.8, group = 1) +
  geom_point(aes(y = n_gene, color = "Total DEGs"), size = 2) +
  
  scale_color_manual(values = c(
    "Higher expressed" = "#B4251A",
    "Lower expressed" = "#163A7D",
    "Total DEGs" = "black"
  )) +
  labs(
    x = "Start of the season boundary",
    y = "Number of significantly expressed genes",
    color = ""
  ) +
  theme(axis.text=element_text(size=15),
      axis.title=element_text(size=15),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom",
      panel.grid.major =  element_line(colour = "grey70", size = 0.1), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black", linewidth =0.6),
      legend.text=element_text(size=12),
      legend.key = element_rect(fill = "white"))

png("./Season/Season_new_date/plots/Figure1_linear.png", width = 1400, 
    height = 1000,
    res = 150)
print(figure1)
dev.off()


volcano <- readRDS("./Season/Season_new_date/plots/Volcano_seasons_15Nov.rds")
#print(volcano)

## combine 2 plots
combined_fig1 <- ggarrange(figure1, volcano,
                           ncol = 1,
                           nrow = 2,
                           heights = c(2.25, 2.5),
                           labels = c("A", "B"),
                           label.x = 1,      # Horizontal position (0 = left, 1 = right)
                           label.y = 1,      # Vertical position (0 = bottom, 1 = top)
                           hjust = 2.5,     # Fine-tune horizontal justification
                           vjust = 2.5) 



png("./Season/Season_new_date/plots/Figure1_linear_volcano.png", 
    width = 1200, 
    height = 1700,
    res = 150)
print(combined_fig1)
dev.off()
