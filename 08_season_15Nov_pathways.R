# this script will perform GSEA analysis for genes higher in winter-spring and lower in winter-spring 
library(data.table)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)

setwd("/Users/tatiana/Work/RP2/ATLANTIS")

set.seed(123)
de.results <- read.csv("./Season/Season_new_date/DE_genes_15Nov_winter_spring.csv")

# convert Hensemble gene IDs to more stable ENTREZ IDs
entrez <- bitr(de.results$Gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
#  13.76% of input gene IDs are fail to map...

# Convert ENTREZ to text, so that R would process it correctly
entrez$ENTREZID <- as.character(entrez$ENTREZID)

# Add column with ENTREZ IDs to the table, remove rows where ENTREZ was missing
de.results <- merge(de.results, entrez, by.x = 'Gene', by.y = 'ENSEMBL') 

#### GSEA with go ######
de.results <- de.results %>%
  mutate(`-10logPval` = ifelse(logFC>0, 1,-1)*-log10(PValue)) %>%
  arrange(desc (`-10logPval`))

up.genes <- de.results %>%
  filter(logFC > 0) %>%
  arrange(desc(`-10logPval`))

up.list <- up.genes$`-10logPval`
names(up.list) <- up.genes$ENTREZID

# gsea_GO_up <- clusterProfiler::gseGO(up.list,
#                                      ont = "ALL",
#                                      keyType = "ENTREZID",
#                                      OrgDb = 'org.Hs.eg.db',
#                                      pvalueCutoff = 0.05,
#                                      minGSSize = 10,
#                                      maxGSSize = 10000,
#                                      scoreType = "pos")
# 
# gsea_GO_up_simp <- clusterProfiler::simplify(gsea_GO_up, cutoff=0.7, by="p.adjust", select_fun=min)
# save(gsea_GO_up_simp, file = "./Season/Season_new_date/GSEA_GO_up_simpl.RData")
load("./Season/Season_new_date/GSEA_GO_up_simpl.RData")
gsea_GO_up_df <- gsea_GO_up_simp@result %>%
  arrange(ONTOLOGY, setSize) %>%
  dplyr::mutate_at(c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"), ~signif(., 3))
write.table(gsea_GO_up_df, "./Season/Season_new_date/GSEA_GO_up_simpl.csv", row.names = F, quote = F, sep = '\t')

down.genes <- de.results %>%
  filter(logFC < 0) %>%
  arrange(desc(`-10logPval`))

down.list <- down.genes$`-10logPval`
names(down.list) <- down.genes$ENTREZID

# gsea_GO_down <- clusterProfiler::gseGO(down.list,
#                                      ont = "ALL",
#                                      keyType = "ENTREZID",
#                                      OrgDb = 'org.Hs.eg.db',
#                                      pvalueCutoff = 0.05,
#                                      minGSSize = 10,
#                                      maxGSSize = 10000,
#                                      scoreType = "neg")
# 
# gsea_GO_down_simp <- clusterProfiler::simplify(gsea_GO_down, cutoff=0.7, by="p.adjust", select_fun=min)
# 
# save(gsea_GO_down_simp, file = "./Season/Season_new_date/GSEA_GO_down_simpl.RData")
load("./Season/Season_new_date/GSEA_GO_down_simpl.RData")

gsea_GO_down_df <- gsea_GO_down_simp@result %>%
  arrange(ONTOLOGY, setSize) %>%
  dplyr::mutate_at(c("enrichmentScore", "NES", "pvalue", "p.adjust", "qvalue"), ~signif(., 3))
write.table(gsea_GO_down_df, "./Season/Season_new_date/GSEA_GO_down_simpl.csv", row.names = F, quote = F, sep = "\t")
