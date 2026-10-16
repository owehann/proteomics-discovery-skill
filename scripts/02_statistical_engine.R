#!/usr/bin/env Rscript
# ==============================================================================
# Stage 2 (R Primary Engine): Bioconductor Statistical Engine for Proteomics
# ==============================================================================
# Supports:
#   1. limpa (Li, Cobbold & Smyth 2025): Detection Probability Curve (dpc / dpcDE)
#   2. DEqMS (Zhu et al. 2020): Peptide/PSM-count weighted empirical Bayes
#   3. limma (Ritchie et al. 2015): Trended/robust eBayes + duplicateCorrelation
#   4. TOST Equivalence Testing: Two One-Sided Tests within [-delta_eq, +delta_eq]
# ==============================================================================

suppressPackageStartupMessages({
  library(methods)
  library(stats)
  library(utils)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript 02_statistical_engine.R <input_csv> <design_csv> <output_csv> [delta_eq_log2]\n")
  quit(status = 0)
}

input_csv  <- args[1]
design_csv <- args[2]
output_csv <- args[3]
delta_eq   <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.321928)

dat <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
des <- read.csv(design_csv, stringsAsFactors = FALSE, check.names = FALSE)

sample_cols <- des$sample
expr_mat <- as.matrix(dat[, sample_cols, drop = FALSE])
rownames(expr_mat) <- dat$composite_id

group <- factor(des$group, levels = unique(des$group))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

has_limpa <- requireNamespace("limpa", quietly = TRUE)
has_limma <- requireNamespace("limma", quietly = TRUE)
has_deqms <- requireNamespace("DEqMS", quietly = TRUE)

if (!has_limma) {
  stop("Bioconductor package 'limma' is required for 02_statistical_engine.R. Use 02_statistical_engine_py.py as fallback.")
}

contrast_str <- paste0(levels(group)[2], " - ", levels(group)[1])
cont_mat <- limma::makeContrasts(contrasts = contrast_str, levels = design)

if (has_limpa) {
  dpc_fit <- tryCatch({
    dpc_obj <- limpa::dpc(expr_mat)
    fit <- limpa::dpcDE(expr_mat, design, dpc = dpc_obj)
    fit2 <- limma::contrasts.fit(fit, cont_mat)
    limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  }, error = function(e) NULL)
} else {
  dpc_fit <- NULL
}

if (is.null(dpc_fit)) {
  fit <- limma::lmFit(expr_mat, design)
  fit2 <- limma::contrasts.fit(fit, cont_mat)
  if (has_deqms && ("std_unique_peptides" %in% colnames(dat))) {
    fit2$count <- pmax(as.integer(dat$std_unique_peptides), 1L)
    fit2 <- tryCatch(
      DEqMS::spectraCounteBayes(fit2),
      error = function(e) limma::eBayes(fit2, trend = TRUE, robust = TRUE)
    )
  } else {
    fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
  }
} else {
  fit2 <- dpc_fit
}

tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
se_vec <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1]
df_total <- fit2$df.total

# Compute Two One-Sided Tests (TOST) for statistical equivalence within [-delta_eq, +delta_eq]
t_lower <- (tt$logFC - (-delta_eq)) / se_vec
t_upper <- (delta_eq - tt$logFC) / se_vec
p_lower <- stats::pt(t_lower, df = df_total, lower.tail = FALSE)
p_upper <- stats::pt(t_upper, df = df_total, lower.tail = FALSE)
tt$tost_p_value <- pmax(p_lower, p_upper)
tt$tost_fdr <- stats::p.adjust(tt$tost_p_value, method = "BH")

out_df <- cbind(
  dat[, intersect(c("composite_id", "std_accession", "std_gene", "std_unique_peptides", "peptide_tier"), colnames(dat)), drop = FALSE],
  tt
)
write.csv(out_df, output_csv, row.names = FALSE)
cat(sprintf("[Stage 2 R Engine] Wrote moderated DE + TOST results (%d features) to %s\n", nrow(out_df), output_csv))
