run_edger_pipeline <- function(expression_data, design) {
  DGEL <- edgeR::DGEList(expression_data)
  keep <- edgeR::filterByExpr(DGEL, design)
  DGEL <- DGEL[keep, , keep.lib.sizes = FALSE]
  DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
  DGEL <- edgeR::estimateDisp(DGEL, design)
  fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) ## reproduce previous version of edgeR
  list(DGEL = DGEL, fit = fit)
}

get_de_results <- function(fit, DGEL, coef = NULL, contrast = NULL, gene_map) {
  qlf <- if (!is.null(contrast)) {
    edgeR::glmQLFTest(fit, contrast = contrast)
  } else {
    edgeR::glmQLFTest(fit, coef = coef)
  }
  print(summary(decideTests(qlf)))

  topTags(qlf, n = nrow(DGEL))$table %>%
    tibble::rownames_to_column("Gene") %>%
    left_join(gene_map, by = c("Gene" = "ensembl_gene_id"))
}

get_gene_names <- function(DGEL) {
  ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  biomaRt::getBM(
    attributes = c("ensembl_gene_id", "external_gene_name"),
    filters = "ensembl_gene_id",
    values = rownames(DGEL$counts),
    mart = ensembl
  )
}