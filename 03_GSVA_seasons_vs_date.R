# This code will calculate GSVA scores for 2 seasonal geneset and plot those scores against the date of the sample
library(dplyr)
library(GSVA)
library(lubridate)
library(ggplot2)
library(edgeR)
library(ggpubr)

setwd("~/Work/RP2/ATLANTIS")

master.Table <- read.csv("./Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv")

# raw expression data
expression.data <- read.csv("./Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv", header =TRUE) %>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))


# Create an edgeR object, filter low expressed genes, normalize
design <- model.matrix(~new_seasons + age + gender + smoking.status + asthma.status, data = master.Table)

DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL, design = design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
ATLANTIS_logcpm <- cpm(DGEL,normalized.lib.sizes=TRUE, log=TRUE)

# upload ATLANTIS season genes and calculate gsva score
DE_ATLANTIS_seasons <- read.csv('./Season/Season_new_date/DE_genes_15Nov_winter_spring.csv')

geneSets <- list(Up_genes = DE_ATLANTIS_seasons[DE_ATLANTIS_seasons$logFC>0 & DE_ATLANTIS_seasons$FDR<0.05,]$Gene,
                 Down_genes = DE_ATLANTIS_seasons[DE_ATLANTIS_seasons$logFC<0 & DE_ATLANTIS_seasons$FDR<0.05,]$Gene)

param.object <- gsvaParam(ATLANTIS_logcpm, geneSets = geneSets)
gsva_ATLANTIS <- gsva(param.object)

gsva_ATLANTIS_trans <- as.data.frame(t(gsva_ATLANTIS))%>%
  tibble::rownames_to_column('Sample')%>%
  left_join(master.Table[,c('GenomeScan_ID','VISIT.DATE', 'new_seasons', 'asthma.status')], by = c('Sample' = 'GenomeScan_ID')) %>%
  mutate(VISIT.DATE = as.Date(VISIT.DATE))

gsva_ATLANTIS_trans$new_seasons <- factor(gsva_ATLANTIS_trans$new_seasons, levels = c("winter_spring", "summer_autumn"))

# plot 

Up <- ggplot(gsva_ATLANTIS_trans, aes(x = VISIT.DATE, y = Up_genes))+
  geom_point(aes(color = new_seasons, shape = asthma.status))+
  geom_smooth(colour = "black")+
  scale_colour_manual(values=c("winter_spring"= "#BFCCB5","summer_autumn"="#FFA559"),
                      aesthetics = "colour",
                      name='',
                      labels = c('winter-spring', 'summer-autumn')) +
  scale_shape_manual(values=c("A"= 1,"H"= 2 ),
                     aesthetics = "shape",
                     name='',
                     labels = c('Asthmatic subjects', 'Healthy subjects')) +
  ylab('Higher expressed genes in winter-spring') +
  xlab('Date') +
  theme_classic() +
  theme(axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) +
  expand_limits(x = as.Date("2016-03-10"))


Down <- ggplot(gsva_ATLANTIS_trans, aes(x = VISIT.DATE, y = Down_genes)) +
  geom_point(aes(color = new_seasons, shape = asthma.status))+
  geom_smooth(colour = "black") +
  scale_colour_manual(values=c("winter_spring"= "#BFCCB5","summer_autumn"="#FFA559"),
                      aesthetics = "colour",
                      name='',
                      labels = c('winter-spring', 'summer-autumn')) +
  scale_shape_manual(values=c("A"= 1,"H"= 2 ),
                     aesthetics = "shape",
                     name='',
                     labels = c('Asthmatic subjects', 'Healthy subjects')) +
  ylab("Lower expressed genes in winter-spring") +
  xlab('Date') +
  theme_classic() +
  theme(axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10)) + 
  expand_limits(x = as.Date("2016-03-10"))

all <- ggarrange(Up + rremove("xlab"), Down,
                 common.legend = TRUE,
                 legend = "right",
                 nrow = 2)
annotated_figure <- annotate_figure(all,
                                    left = text_grob("GSVA value",  rot = 90, face = "bold"))
png("./Season/Season_new_date/plots/Date_vs_gsva_up_down.png", width=1400, height=1300,res = 150)
print(annotated_figure)
dev.off()

saveRDS(all, "./Season/Season_new_date/plots/Date_vs_gsva_up_down.RDS")
