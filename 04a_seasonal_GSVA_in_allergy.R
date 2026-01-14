# This script will calculate GSVA score for seasonal gene sets and compare it between allergic and  non-allergic individuals 
library(dplyr)
library(GSVA)
library(lubridate)
library(ggplot2)
library(rstatix)
library(ggpubr)

setwd("~/Work/RP2/ATLANTIS/")

master.Table <- read.csv("./Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv")

# Add allergy data 

clinical_table <-  read.csv('./atlantis_patient_data.csv', header =TRUE, na.strings=c("","NA"))%>%
  #dplyr::select(-c(X))%>%
  dplyr::select(c('PT','PHADRES')) 
clinical_table <- clinical_table[!duplicated(clinical_table$PT), ]

master.Table.ATLANTIS <- master.Table%>%
  left_join(clinical_table, by = ('PT' = 'PT'))

master.Table.short <- master.Table.ATLANTIS %>%
  filter(!is.na(PHADRES)) %>%
  mutate(PHADRES = as.factor(PHADRES))

master.Table.short%>%
  group_by(PHADRES,asthma.status) %>%
  summarise(n = n())

# plot dates of sampling
ggplot(master.Table.short, aes (x =as.Date(VISIT.DATE) , y = PHADRES))+
  geom_point(aes(color = PHADRES))+
  theme_minimal()

# raw data
expression.data <- read.csv('./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv', header =TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table.short$GenomeScan_ID))

# normalized data
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
ATLANTIS_logcpm <- edgeR::cpm(DGEL,normalized.lib.sizes=TRUE, log=TRUE)

# upload ATLANTIS season genes and calculate gsva score
DE_ATLANTIS_seasons <- read.csv('./Season/Season_new_date/DE_genes_15Nov_winter_spring.csv')

geneSets <- list(Up_genes = DE_ATLANTIS_seasons[DE_ATLANTIS_seasons$logFC>0 & DE_ATLANTIS_seasons$FDR<0.05,]$Gene,
                 Down_genes = DE_ATLANTIS_seasons[DE_ATLANTIS_seasons$logFC<0 & DE_ATLANTIS_seasons$FDR<0.05,]$Gene)

param.object <- gsvaParam(ATLANTIS_logcpm, geneSets = geneSets)
gsva_ATLANTIS <- gsva(param.object)

gsva_ATLANTIS_trans <- as.data.frame(t(gsva_ATLANTIS))%>%
  tibble::rownames_to_column('Sample')%>%
  left_join(master.Table.short[,c('GenomeScan_ID','VISIT.DATE', 'new_seasons', 'asthma.status', 'PHADRES')],
            by = c('Sample' = 'GenomeScan_ID')) %>%
  mutate(VISIT.DATE = as.Date(VISIT.DATE))

########### plot ################
df_ATLANTIS_GSVA <- gsva_ATLANTIS_trans%>%
  select(c('Sample', 'Up_genes','Down_genes'))%>%
  tidyr::gather(direction, GSVA_value, -Sample)%>%
  left_join(gsva_ATLANTIS_trans%>%
              select(c('Sample','new_seasons', 'PHADRES')))%>%
  mutate(PHADRES = if_else(PHADRES == 'Negative', 'Non-allergic subjects', 'Allergic subjects'))

Up_genes_plot <- df_ATLANTIS_GSVA%>%
  filter(direction == 'Up_genes')

Up_genes_plot$new_seasons <- factor(
  Up_genes_plot$new_seasons,
  levels = c("summer_autumn", "winter_spring"))

stat.test.up <- Up_genes_plot%>%
  group_by(PHADRES)%>%
  rstatix::wilcox_test(GSVA_value~new_seasons)%>%
  add_y_position(step.increase = 0.08) %>%
  add_x_position(x = "PHADRES", dodge = 0) %>%
  add_significance("p",
                   cutpoints = c(0, 0.01, 0.05, 1),
                   symbols = c("<0.01", "<0.05", "ns")) %>%
  mutate(p_label = I(sapply(p, function(x) {
    sci <- formatC(x, format = "e", digits = 2)   # e.g. "8.90e-07"
    parts <- strsplit(sci, "e")[[1]]
    mantissa <- parts[1]
    exponent <- as.integer(parts[2])
    bquote("p = " ~ .(mantissa) %*% 10^.(exponent))
  })))

up_allergy <- ggplot(Up_genes_plot, aes(x = PHADRES, y = GSVA_value))+
  geom_boxplot(aes(fill = new_seasons),outlier.shape=NA)+
  geom_point(aes(fill = new_seasons), shape = 21, position=position_jitterdodge(0.2), alpha= 0.5)+
  ylab(expression('Higher expressed genes in winter-spring'))+
  xlab(expression(bold('Allergy status')))+
  #ggtitle("Seasonal changes for allergic and non-allergic groups")+
  scale_fill_manual(values=c("#FFA559","#BFCCB5"), 
                    aesthetics = "fill", 
                    name='',
                    labels = c('summer-autumn','winter-spring')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size =0.6),
        axis.title.x=element_blank(),
        axis.text = element_text(face = "plain", size = 10, colour = 'black'),
        #axis.title = element_text(face = "plain", size = 10, colour = 'black'),
        plot.title = element_text(color="black", size = 14, face = "plain", hjust = 0.5))+
  stat_pvalue_manual(stat.test.up, label = "p_label", xmin = 'xmin', xmax ='xmax',size = 3,
                     remove.bracket = TRUE, parse = TRUE) +
  expand_limits(y = max(Up_genes_plot$GSVA_value) * 1.1)


save(up_allergy, file = "./Season/Season_new_date/plots/up_allergy_seasons_GSVA.Rdata")


Down_genes_plot <- df_ATLANTIS_GSVA%>%
  filter(direction == 'Down_genes')

Down_genes_plot$new_seasons <- factor(
  Down_genes_plot$new_seasons,
  levels = c("summer_autumn", "winter_spring"))

stat.test.down <- Down_genes_plot%>%
  group_by(PHADRES)%>%
  rstatix::wilcox_test(GSVA_value~new_seasons)%>%
  add_y_position(step.increase = 0.08) %>%
  add_x_position(x = "PHADRES", dodge = 0) %>%
  add_significance("p",
                   cutpoints = c(0, 0.01, 0.05, 1),
                   symbols = c("<0.01", "<0.05", "ns")) %>%
  mutate(p_label = I(sapply(p, function(x) {
    sci <- formatC(x, format = "e", digits = 2)   # e.g. "8.90e-07"
    parts <- strsplit(sci, "e")[[1]]
    mantissa <- parts[1]
    exponent <- as.integer(parts[2])
    bquote("p = " ~ .(mantissa) %*% 10^.(exponent))
  })))

down_allergy <- ggplot(Down_genes_plot, aes(x = PHADRES, y = GSVA_value))+
  geom_boxplot(aes(fill = new_seasons),outlier.shape=NA)+
  geom_point(aes(fill = new_seasons), shape = 21, position=position_jitterdodge(0.2), alpha= 0.5)+
  ylab(expression('Lower expressed genes in winter-spring'))+
  xlab(expression(bold('Allergy status')))+
  scale_fill_manual(values=c("#FFA559","#BFCCB5"), 
                    aesthetics = "fill", 
                    name='',
                    labels = c('summer-autumn','winter-spring')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size =0.6),
        axis.title.x=element_blank(),
        axis.text = element_text(face = "plain", size = 10, colour = 'black'),
        axis.title = element_text(face = "plain", size = 10, colour = 'black'))+
  stat_pvalue_manual(stat.test.down, label = "p_label", xmin = 'xmin', xmax ='xmax',size = 3,
                     remove.bracket = TRUE, parse = TRUE) +
  expand_limits(y = max(Down_genes_plot$GSVA_value) * 1.1)

save(down_allergy, file = "./Season/Season_new_date/plots/down_allergy_seasons_GSVA.Rdata")


all_allergy <- ggarrange(up_allergy,down_allergy,
                         #labels = c("A", "B"),
                         common.legend = TRUE,
                         legend = "bottom",
                         nrow = 2)
save(all_allergy, file = "./Season/Season_new_date/plots/Allergy_seasons_GSVA.Rdata")


