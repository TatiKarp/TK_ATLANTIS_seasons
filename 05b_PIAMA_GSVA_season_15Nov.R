# This script will a) calculate GSVA score for seasonal gene sets 
# b) compare the score values between the seasons
# c) plot it against the data of the sample collection
library(dplyr)
library(haven)
library(lubridate)
library(ggplot2)
library(ggpubr)
library(tibble)
library(GSVA)
library(matrixStats)
library(ggpubr)
library(rstatix)
library(config)

# location 
conf <- config::get()

# load data
PIAMA_gsva_raw <- read.csv(file.path(conf$data_path, "Season/replication_PIAMA/PIAMA_gsva_scores.csv"), row.names = 'X') %>%
  tibble::rownames_to_column('Sample')

PIAMA_dates <- read_sav(file.path(conf$data_path, "Season/replication_PIAMA/datesage16.sav")) %>%
  mutate(ID = as.character(ID))

pheno_PIAMA <- read.csv(file.path(conf$data_path, "Season/replication_PIAMA/piama_rnaseq_subjects.csv"), sep = ';') %>%
  mutate(ID_numb = substring(subject_ID, 1,5))

### update IDs 
piama_ts_newid<- function(sampleid)
  ## transform old_sampleid to new_sampleid
{
  key<- read.csv(file.path(conf$data_path, "Season/replication_PIAMA/Keyid.csv"))
  nindex<-seq(1,length(sampleid))
  index<- match(sampleid,key[,1])
  index.c<- cbind(index,nindex); 
  index.c2<- na.omit(index.c)
  newid<- sampleid
  newid[index.c2[,2]]<- key[index.c2[,1],2] 
  return (newid)
}

PIAMA_gsva <- PIAMA_gsva_raw %>%
  mutate(RNA_IDs_numb = substring(Sample,2,6))
PIAMA_gsva$RNA_IDs_new <- piama_ts_newid(PIAMA_gsva$RNA_IDs_numb)

#### add dates to gsva
PIAMA_gsva <- PIAMA_gsva%>%
  left_join(PIAMA_dates%>%
              dplyr::select(c(ID, pdate)), by = c('RNA_IDs_new'='ID'))

#### add asthma to gsva
PIAMA_gsva <- PIAMA_gsva%>%
  left_join(pheno_PIAMA%>%
              dplyr::select(c(ID_numb, asthma.status)), by = c('RNA_IDs_new'='ID_numb'))%>%
  mutate(asthma.status = as.factor(asthma.status)) %>%
  mutate(date_yday = yday(pdate)) %>%
  mutate(new_seasons = as.factor(if_else(((date_yday >= 135) & (date_yday < 319)), "summer_autumn", "winter_spring")),
         week_day = weekdays(pdate))  

#### upload ATLANTIS season genes
DE_ATLANTIS_seasons<-read.csv(file.path(conf$data_path, "Season/Season_new_date/DE_genes_15Nov_winter_spring.csv"))


#### determine gene sets (up and down)
geneSets <- list(Up_genes = DE_ATLANTIS_seasons[DE_ATLANTIS_seasons$logFC>0 & DE_ATLANTIS_seasons$FDR<0.05,]$Gene,
                 Down_genes = DE_ATLANTIS_seasons[DE_ATLANTIS_seasons$logFC<0 & DE_ATLANTIS_seasons$FDR<0.05,]$Gene)


master.Table <- read.csv(file.path(conf$data_path, "Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv"))

expression.data.ATLANTIS <- read.csv(file.path(conf$data_path, "Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header =TRUE) %>%
  tibble::column_to_rownames("Gene") %>%
  dplyr::select(c(master.Table$GenomeScan_ID))%>%
  as.matrix()


DGEL_ATLANTIS <- edgeR::DGEList(expression.data.ATLANTIS)
keep <- edgeR::filterByExpr(DGEL_ATLANTIS, design=NULL) 
DGEL_ATLANTIS <- DGEL_ATLANTIS[keep, , keep.lib.sizes=FALSE]
DGEL_ATLANTIS <- edgeR::calcNormFactors(DGEL_ATLANTIS, method = "TMM")
logcpm_ATLANTIS <- cpm(DGEL_ATLANTIS, normalized.lib.sizes = TRUE, log = TRUE)%>%
  as.matrix()

###### GSVA ATLANTIS #######
param.object <- gsvaParam(logcpm_ATLANTIS, geneSets = geneSets)
gsva_ATLANTIS <- gsva(param.object)

gsva_ATLANTIS_trans<-as.data.frame(t(gsva_ATLANTIS))%>%
  rownames_to_column('Sample')%>%
  tidyr::gather(direction, GSVA_value, -Sample)%>%
  mutate(Study_name =  'ATLANTIS')%>%
  left_join(master.Table[,c('GenomeScan_ID','new_seasons', 'VISIT.DATE')], by = c('Sample' = 'GenomeScan_ID'))

######### all PIAMA ##############
gsva_PIAMA_trans<-PIAMA_gsva_raw%>%
  tidyr::gather(direction, GSVA_value, -Sample)%>%
  mutate(Study_name =  'PIAMA')%>%
  left_join(PIAMA_gsva[,c('Sample','new_seasons', 'pdate')], by = c('Sample' = 'Sample'))%>%
  mutate(pdate = as.character(pdate))%>%
  rename(Date = pdate)%>%
  tidyr::drop_na()


GSVA_ATLANTIS_PIAMA <-bind_rows(gsva_ATLANTIS_trans, gsva_PIAMA_trans)

################## PLOT ##################

Up_genes_plot <- GSVA_ATLANTIS_PIAMA%>%
  filter(direction == 'Up_genes')%>%
  arrange(new_seasons)
Up_genes_plot$new_seasons <- factor(
  Up_genes_plot$new_seasons, 
  levels = c("summer_autumn", "winter_spring"))


stat.test.up <- Up_genes_plot%>%
  group_by(Study_name)%>%
  rstatix::wilcox_test(GSVA_value~new_seasons)%>%
  add_y_position(step.increase = 0.08) %>%
  add_x_position(x = "Study_name", dodge = 0) %>%
  add_significance("p",
                   cutpoints = c(0, 0.01, 0.05, 1),
                   symbols = c("<0.01", "<0.05", "ns")) %>%
  mutate(p.signif = if_else(p.signif == "ns", paste0("p=", as.character(round(p,2))), p.signif)) %>%
  mutate(p_label = I(sapply(p, function(x) {
    sci <- formatC(x, format = "e", digits = 2)   # e.g. "8.90e-07"
    parts <- strsplit(sci, "e")[[1]]
    mantissa <- parts[1]
    exponent <- as.integer(parts[2])
    bquote("p = " ~ .(mantissa) %*% 10^.(exponent))
  })))
  

up <- ggplot(Up_genes_plot, aes(x = Study_name, y = GSVA_value))+
  geom_boxplot(aes(fill = new_seasons),outlier.shape=NA)+
  geom_point(aes(fill = new_seasons), shape = 21, position=position_jitterdodge(0.2), alpha = 0.5)+
  #stat_summary(fun="median", geom = 'text', aes(label=round(..y.., digits=4), group= asthma.status),color = 'black',
  #             position=position_dodge(width = 0.8), vjust=-0.5, size = 3)+
  #stat_compare_means(method = "wilcox.test", paired = FALSE, label = 'p.format', label.y = 1)+
  ylab(expression("Higher expressed genes in winter-spring"))+
  xlab(' ')+
  scale_fill_manual(values=c("#FFA559","#BFCCB5"), aesthetics = "fill", name='',labels = c('summer-autumn','winter-spring'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size =0.6),
        axis.text = element_text(face = "plain", size = 10, colour = 'black'),
        axis.title = element_text(face = "plain", size = 10, colour = 'black'))+
  stat_pvalue_manual(stat.test.up, label = "p_label", xmin = 'xmin', xmax ='xmax',size = 3,
                     remove.bracket = TRUE, parse = TRUE) +
  expand_limits(y = max(Up_genes_plot$GSVA_value) * 1.1)


save(up, file = file.path(conf$data_path, "Season/Season_new_date/plots/up_PIAMA_ATLANTIS_seasons_GSVA.Rdata"))

### plot down genes:

Down_genes_plot <- GSVA_ATLANTIS_PIAMA%>%
  filter(direction == 'Down_genes')
Down_genes_plot$new_seasons <- factor(
  Down_genes_plot$new_seasons, 
  levels = c("summer_autumn", "winter_spring"))

stat.test.down <- Down_genes_plot%>%
  group_by(Study_name)%>%
  rstatix::wilcox_test(GSVA_value~new_seasons)%>%
  add_y_position(step.increase = 0.08) %>%
  add_x_position(x = "Study_name", dodge = 0) %>%
  add_significance("p",
                   cutpoints = c(0, 0.01, 0.05, 1),
                   symbols = c("<0.01", "<0.05", "ns")) %>%
  mutate(p.signif = if_else(p.signif == "ns", paste0("p=", as.character(round(p,2))), p.signif)) %>%
  mutate(p_label = I(sapply(p, function(x) {
    sci <- formatC(x, format = "e", digits = 2)   # e.g. "8.90e-07"
    parts <- strsplit(sci, "e")[[1]]
    mantissa <- parts[1]
    exponent <- as.integer(parts[2])
    bquote("p = " ~ .(mantissa) %*% 10^.(exponent))
  })))


down <- ggplot(Down_genes_plot, aes(x = Study_name, y = GSVA_value))+
  geom_boxplot(aes(fill = new_seasons),outlier.shape=NA)+
  geom_point(aes(fill = new_seasons), shape = 21, position=position_jitterdodge(0.2), alpha= 0.5)+
  #stat_summary(fun="median", geom = 'text', aes(label=round(..y.., digits=4), group= asthma.status),color = 'black',
  #             position=position_dodge(width = 0.8), vjust=-0.5, size = 3)+
  #stat_compare_means(method = "wilcox.test", paired = FALSE, label = 'p.format', label.y = 1)+
  #ggtitle('Downregulated genes')+
  ylab(expression("Lower expressed genes in winter-spring"))+
  xlab(expression("Study name"))+
  scale_fill_manual(values=c("#FFA559","#BFCCB5"), aesthetics = "fill", name='',labels = c("summer-autumn", "winter-spring"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size =0.6),
        axis.text = element_text(face = "plain", size = 10, colour = 'black'),
        axis.title = element_text(face = "plain", size = 10, colour = 'black')) +
  stat_pvalue_manual(stat.test.down, label = "p_label", xmin = 'xmin', xmax ='xmax',size = 3,
                     remove.bracket = TRUE, parse = TRUE) +
  expand_limits(y = max(Down_genes_plot$GSVA_value) * 1.1)


save(down, file = file.path(conf$data_path, "Season/Season_new_date/plots/down_PIAMA_ATLANTIS_seasons_GSVA.Rdata"))

up_down <- ggarrange(up,down,
                     labels = c("A", "B"),
                     common.legend = TRUE,
                     legend = "bottom",
                     nrow = 2)
annotated_figure <- annotate_figure(up_down,
                                    left = text_grob("GSVA value",  rot = 90, face = "bold"))
png(file.path(conf$data_path, "Season/Season_new_date/plots/GSVA_piama_ATLANTIS_all.png"), width = 1100, height = 1200,res = 150)
print(annotated_figure)
dev.off()


#####  GSVA over the year in PIAMA #####
PIAMA_gsva <- PIAMA_gsva %>%
  filter(!is.na(pdate))

up <- ggplot(PIAMA_gsva, aes (x = pdate, y = Up_genes))+
  geom_point(aes(color = new_seasons, shape = asthma.status))+
  geom_smooth(colour = "black")+
  scale_colour_manual(values=c("winter_spring"= "#BFCCB5","summer_autumn"="#FFA559"), 
                      aesthetics = "colour",
                      name='',
                      labels = c('summer-autumn', 'winter-spring')) +
  scale_shape_manual(values=c("asthma"= 1,"no asthma"= 2 ),
                     aesthetics = "shape",
                     name='',
                     labels = c('Asthmatic subjects', 'Healthy subjects')) +
  ylab("Higher expressed genes in winter-spring") +
  xlab("Date") + 
  scale_x_date(date_labels = "%Y-%m") +
  theme_classic() +
  theme(axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))



down <- ggplot(PIAMA_gsva, aes (x = pdate, y = Down_genes))+
  geom_point(aes(color = new_seasons, shape = asthma.status))+
  geom_smooth(colour = "black")+
  scale_colour_manual(values=c("winter_spring"= "#BFCCB5","summer_autumn"="#FFA559"), 
                      aesthetics = "colour",
                      name='',
                      labels = c('summer-autumn', 'winter-spring')) +
  scale_shape_manual(values=c("asthma"= 1,"no asthma"= 2 ),
                     aesthetics = "shape",
                     name='',
                     labels = c('Asthmatic subjects', 'Healthy subjects')) +
  ylab("Lower expressed genes in winter-spring") +
  xlab("Date") +
  scale_x_date(date_labels = "%Y-%m") +
  theme_classic() +
  theme(axis.text = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_text(size=10))

all <- ggarrange(up + rremove("xlab"), down,
                 common.legend = TRUE,
                 legend = "bottom",
                 nrow = 2)

save(all, file = file.path(conf$data_path, "Season/Season_new_date/plots/Date_vs_gsva_up_down_PIAMA.Rdata"))

annotated_figure <- annotate_figure(all,
                                    left = text_grob("GSVA score",  rot = 90, face = "bold"))
png(file.path(conf$data_path, "Season/Season_new_date/plots/PIAMA_GSVA_vs_Date.png"), width = 1500, height = 1400, res = 150)
print(annotated_figure)
dev.off()
