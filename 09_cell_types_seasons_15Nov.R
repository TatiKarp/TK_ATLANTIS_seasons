# This script will  a) perform cell type deconvolution 
# b) compare cell types proportion between the 2 seasons
library(dplyr)
library(tibble)
library(ggplot2)
library(lubridate)
library(config)

# location 
conf <- config::get()

## upload source Cibersort
source(file.path(conf$data_path, "Deconvolution/Deconvolution_Jos/CIBERSORT.R"), verbose=TRUE)

## upload signature file 
# matrix from Tessa 
C <- read.csv(file.path(conf$data_path, "Season/cell_types/CIBERSORTx_nose_subsampled_matrix_max200cells_ENSG_inferred_phenoclasses.CIBERSORTx_nose_subsampled_matrix_max200cells_ENSG_inferred_refsample.bm.K999.txt"), sep='\t')
rownames(C)<-C$NAME
C <- C%>%
  dplyr::select(-NAME)

## bulk data preparation
master.Table <- read.csv(file.path(conf$data_path, "Season/Season_new_date/ATLANTIS_master_table_seasons_15Nov.csv"))

expression.data <- read.csv(file.path(conf$data_path, "Umi_dedup/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header =TRUE) %>%
  tibble::column_to_rownames("Gene")%>%
  dplyr::select(c(master.Table$GenomeScan_ID))

## cpm normalization
bulk = edgeR::cpm(expression.data)
bulk <- bulk %>%
  as.data.frame() %>%
  tibble::rownames_to_column('Gene')
# left_join(gene.data %>%
#             dplyr::select(c(hgnc_symbol, ensembl_gene_id)), by = c('Gene'='ensembl_gene_id'))%>%
# dplyr::select(-c(Gene))

## overlap genes from signature and from bulk data
keep = intersect(rownames(C), bulk$Gene) ## keep 2034 out of the 2035

bulk <- bulk %>%
  filter(Gene %in% keep) %>%
  column_to_rownames("Gene")

Ref = C[keep,]

##############################################################
##################### run Cibersort ##########################

xf <- file.path(conf$data_path, "Season/cell_types/reference.tsv")
Ref %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rownames")%>%
  readr::write_tsv(
    file = xf
  )

yf <-  file.path(conf$data_path, "Season/cell_types/mixture.tsv")
bulk %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rownames") %>%
  readr::write_tsv(
    file = yf
  )

RESULTS <- CIBERSORT(sig_matrix = xf, mixture_file = yf, QN = FALSE, perm=100)
res.cibersort <- t(RESULTS[,1:(ncol(RESULTS)-3)]) %>%
  as.data.frame()

write.csv(res.cibersort, file.path(conf$data_path, "Season/cell_types/CIBERSORT.proportions.csv"), row.names = T, quote = F)

##############################################################
## check edgeR normalization! for the CIBERSORT

CIBESORT_ATL <- read.csv(file.path(conf$data_path, "Season/cell_types/CIBERSORT.proportions.csv")) %>%
  column_to_rownames("X") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample")
master.Table <- master.Table %>%
  mutate(Date = as.Date(VISIT.DATE, format = "%Y-%m-%d"))%>%
  left_join(CIBESORT_ATL, by = c("GenomeScan_ID" = "Sample"))


## plot date vs cell type proportions 

ggplot(master.Table, aes (x = Date, y = Goblet))+
  geom_point()+
  geom_smooth()+
  ylab('Goblet cells')+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(master.Table, aes (x = Date, y = `Multiciliated.lineage`))+
  geom_point()+
  geom_smooth()+
  ylab('Ciliated cells')+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

ggplot(master.Table, aes (x = Date, y = `Basal.resting...Suprabasal`))+
  geom_point()+
  geom_smooth()+
  ylab('Basal.resting cells')+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

## plot boxplots 
library(rstatix)
library(ggpubr)
library(stringr)


CIBESORT_ATL_long <- CIBESORT_ATL %>%
  tidyr::pivot_longer(-Sample,
                      names_to = "Cell_type",
                      values_to = "Proportion") %>%
  left_join(master.Table %>%
              dplyr::select(c(GenomeScan_ID, new_seasons)), by = c("Sample" = "GenomeScan_ID")) %>%
  mutate(Cell_type = case_when (
    Cell_type == "Basal.resting...Suprabasal" ~ "Basal",
    Cell_type == "Dendritic.cells" ~ "Dendritic", 
    Cell_type == "Hillock.like" ~ "Hillock",
    Cell_type == "Multiciliated.lineage" ~ "Ciliated", 
    Cell_type == "T.cell.lineage" ~ "T cells",
    .default = Cell_type))
# mutate(Cell_type = str_replace_all(Cell_type, "\\.", " ")) %>%
# mutate(Cell_type = str_replace_all(Cell_type, "   ", " ")) %>%
# mutate(Cell_type = if_else(Cell_type == "Basal resting Suprabasal", "Basal resting, Suprabasal", Cell_type))

stat.test.CIBER_season <- CIBESORT_ATL_long%>%
  group_by(Cell_type)%>%
  rstatix::wilcox_test(Proportion~new_seasons, detailed = TRUE)%>%
  add_y_position(step.increase = 0.08) %>%
  add_x_position(x = "Cell_type", dodge = 0)%>%
  mutate(p.adj = round(p.adjust(p, method = 'fdr'), digits = 3))

stat.test.sign <- stat.test.CIBER_season%>%
  filter(p < 0.05)

CIBESORT_ATL_long$new_seasons <- factor(CIBESORT_ATL_long$new_seasons, levels = c("summer_autumn", "winter_spring"))

## plot has to be fixed !!!!!
bxp_seasons <- CIBESORT_ATL_long %>%
  ggplot(mapping = aes(x = Cell_type, y = Proportion)) +
  geom_boxplot(aes(fill = new_seasons), outlier.shape = NA)+
  ylab('Proportion')+
  xlab('Cell type')+
  #stat_compare_means(method = "wilcox.test", paired = FALSE, label = 'p.format', label.y = 1)+
  scale_fill_manual(values=c("#FFA559", "#BFCCB5"), aesthetics = "fill", name='', labels = c("summer-autumn", "winter-spring"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black", size =0.6),
        axis.text = element_text(face = "plain", size = 10, colour = 'black'),
        axis.title = element_text(face = "plain", size = 10, colour = 'black'),
        axis.text.x = element_text(angle = 60, hjust = 1))+
  geom_point(aes(fill = new_seasons), shape = 21, position = position_jitterdodge(jitter.width = 0.2), alpha = 0.3) +
  # stat_pvalue_manual(stat.test.CIBER_season, label = "p ={p}", xmin = 'xmin', xmax ='xmax',size = 3,
  #                    remove.bracket = TRUE, y.position = stat.test.CIBER_season$y.position + 0.04) +
  stat_pvalue_manual(stat.test.CIBER_season, label = "p_adj = {p.adj}", xmin = 'xmin', xmax ='xmax',size = 3, remove.bracket = TRUE)
#annotate("text", x = Inf, y = 1, hjust = 1, vjust = 0, label = "wilcox-test, FDR adjusted")

png(file.path(conf$data_path, "Season/Season_new_date/plots/CIBERSORT_seasons_15Nov.png"), 
    res = 150, 
    width = 1200, 
    height = 800)
print(bxp_seasons)
dev.off()

save(bxp_seasons, file = file.path(conf$data_path, "Season/Season_new_date/plots/CIBERSORT_seasons_15Nov.Rdata"))


## get medians, IQR, Q1 and Q3 from the table 
CIBESORT_ATL_long %>%
  group_by(Cell_type, new_seasons) %>%
  summarize(median_value = median(Proportion),
            iqr_value = IQR(Proportion),
            q1_value = quantile(Proportion, 0.25),
            q3_value = quantile(Proportion, 0.75))

