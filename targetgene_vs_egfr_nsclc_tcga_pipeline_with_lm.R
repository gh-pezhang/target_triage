###########################################################################
## script to use TCGA data (LUAD + LUSC) to evaluate candidate targets
## to overcome osimertinib resistance in EGFR+ NSCLC, by following analyses
## 1. expression correlation with EGFR
## 2. expression correlation with bypass genes
## 3. EGFR-driver-mutant vs WT tumor expression: log2FC + wilcoxon p-value
## 4. build linear models between target expression ~ EGFR-mutant subtype
## generate scatter plot, volcano plots and csv tables for the analysis
##
##
## EGFR-driven NSCLC TCGA pipeline
## LUAD + LUSC, MC3 oncogenic EGFR drivers
## + linear models between target expression and EGFR-mutant subtype
## + volcano plot (log2FC vs p-value) for candidate genes
###########################################################################

## 0) User parameters ###########################################
setwd("/mnt/target_triage")

# Path to MC3 MAF file (download separately from GDC / Broad / UCSC)
mc3_maf_file <- "/domino/datasets/local/cohort_comparisons/TCGA/mc3.v0.2.8.PUBLIC.maf"

# Output directories
output_dir <- "./osi_resistance_targets_tcga_results"
plots_dir  <- file.path(output_dir, "plots")
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)

# Target genes triaging (29)
target_genes <- c(
  "IL16", "PTPN6", "IL1A", "PLAT", "LY96", "PTH1R", "ID3", "CD7", 
  "PTEN", "CD47", "GALM", "IL18R1", "RHOH", "CARTPT", "IL2RB", 
  "ADORA2B", "NT5E", "ICOSLG", "MET", "CDKN1C", "TNFRSF12A", 
  "PFKP", "AXL", "HK2", "CELSR2", "DYRK3", "CYP39A1", "LTBR", "HDAC9"
)

# Bypass pathway genes to check co-expression with
bypass_genes <- c(
  "MET","ERBB2","ERBB3","ERBB4",
  "AXL","FGFR1","FGFR2","FGFR3",
  "IGF1R"
)

egfr_gene <- "EGFR"

## 1) Packages ###################################################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("TCGAbiolinks","SummarizedExperiment","DESeq2","maftools")
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
}

cran_pkgs <- c("dplyr","tibble","ggplot2","stringr","readr", "ggrepel")
for (p in cran_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

outdated_packages <- old.packages()
if (!requireNamespace("tidyselect", quietly = TRUE) | 'tidyselect' %in% outdated_packages) {
  install.packages("tidyselect")
}

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(maftools)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(readr)
library(ggrepel)

################ helper functions  ################
map_ensembl_to_hgnc <- function(count_matrix, mart = NULL) {
  # count_matrix: numeric matrix with Ensembl gene IDs as rownames
  # mart: optional biomaRt mart object; if NULL, will create one
  
  if (is.null(mart)) {
    if (!requireNamespace("biomaRt", quietly = TRUE))
      install.packages("biomaRt")
    library(biomaRt)
    
    mart <- useEnsembl(
      biomart = "genes",
      dataset = "hsapiens_gene_ensembl"
    )
  }
  
  # Extract Ensembl IDs and remove version suffixes
  ensembl_ids_raw <- rownames(count_matrix)
  ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids_raw)
  
  # Query BioMart for mapping
  mapping <- getBM(
    filters = "ensembl_gene_id",
    attributes = c(
      "ensembl_gene_id",
      "hgnc_symbol",
      "entrezgene_id"
    ),
    values = ensembl_ids_clean,
    mart = mart
  )
  
  # Remove duplicated Ensembl IDs (keep first occurrence)
  mapping <- mapping[!duplicated(mapping$ensembl_gene_id), ]
  
  # Build lookup vectors
  symbol_lookup <- mapping$hgnc_symbol
  names(symbol_lookup) <- mapping$ensembl_gene_id
  
  entrez_lookup <- mapping$entrezgene_id
  names(entrez_lookup) <- mapping$ensembl_gene_id
  
  # Map Ensembl → HGNC symbol
  mapped_symbols <- symbol_lookup[ensembl_ids_clean]
  
  # Replace empty or NA symbols with original Ensembl ID
  final_symbols <- ifelse(
    is.na(mapped_symbols) | mapped_symbols == "",
    ensembl_ids_clean,
    mapped_symbols
  )
  
  # Apply new rownames to the count matrix
  rownames(count_matrix) <- final_symbols
  
  # Create returned annotation table
  annotation <- data.frame(
    ensembl_id = ensembl_ids_clean,
    hgnc_symbol = mapped_symbols,
    entrez_id = entrez_lookup[ensembl_ids_clean],
    stringsAsFactors = FALSE
  )
  
  return(list(
    counts_mapped = count_matrix,
    annotation = annotation
  ))
}

## 2) Load MC3 MAF and identify EGFR oncogenic/likely-oncogenic ################
## NOTE: MC3 MAF does NOT have tcga_project; we use barcodes only here
## and rely on RNA-seq metadata later to restrict to LUAD/LUSC.

message("Reading MC3 MAF...")
maf_mc3 <- read.maf(mc3_maf_file)

maf_df <- maf_mc3@data

#if (!"ONCOGENIC" %in% colnames(maf_df)) {
#  stop("MC3 MAF does not contain ONCOGENIC annotation column – check that you have MC3+Oncotator version.")
#}
## CLIN_SIG column in maf_df annotates if the mutation is pathogenic

# EGFR oncogenic / likely-oncogenic drivers (all TCGA; cancer-type filtering happens later)
egfr_drivers <- maf_df %>%
  dplyr::filter(
    Hugo_Symbol == egfr_gene &
      CLIN_SIG %in% c("pathogenic","likely_pathogenic", "drug_response")
  ) %>%
  unique()

message("Number of EGFR driver mutations (all TCGA): ", nrow(egfr_drivers))

egfr_driver_patients <- substr(egfr_drivers$Tumor_Sample_Barcode, 1, 12) %>%
  unique()

message("Number of EGFR-driven patients (all TCGA): ", length(egfr_driver_patients))

# Save list of EGFR-driven patients
write_tsv(
  tibble(patient_id = egfr_driver_patients),
  file.path(output_dir, "tcga_egfr_driver_patient_ids.tsv")
)

## 3) Download and prepare RNA-seq for TCGA-LUAD + TCGA-LUSC ###################

projects <- c("TCGA-LUAD","TCGA-LUSC")
se_list <- list()

for (prj in projects) {
  message("Querying RNA-seq for ", prj, " ...")
  GDC_dir <- paste('/domino/datasets/local/cohort_comparisons/TCGA/', prj, sep="")
  
  q <- GDCquery(
    project       = prj,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(q, directory = GDC_dir)
  se_list[[prj]] <- GDCprepare(q, directory = GDC_dir)
}

# Combine SummarizedExperiments
message("Combining LUAD and LUSC SummarizedExperiments...")
# Extract counts and metadata
count_list <- lapply(se_list, function(se) assay(se))
meta_list  <- lapply(se_list, function(se) as.data.frame(colData(se)))

# Align genes across projects
common_genes <- Reduce(intersect, lapply(count_list, rownames))
count_list   <- lapply(count_list, function(mat) mat[common_genes, ])

# Combine counts (genes x samples)
counts_all <- do.call(cbind, count_list)

# Combine metadata using dplyr::bind_rows (unions columns)
meta_all <- bind_rows(meta_list)

# Keep only samples present in counts_all and order rows to match column order
meta_all <- meta_all[match(colnames(counts_all), meta_all$barcode), ]

# Optional: sanity check
stopifnot(all(meta_all$barcode == colnames(counts_all)))

se_all <- SummarizedExperiment(
  assays  = list(counts = counts_all),
  colData = meta_all
)

## 4) Construct EGFR-driver status in expression metadata #######################
## and restrict to LUAD/LUSC via project_id (since MC3 has no project)

# project_id here is from RNA-seq metadata: "TCGA-LUAD" or "TCGA-LUSC"
table(meta_all$project_id)

meta_all$patient_id <- substr(meta_all$barcode, 1, 12)

meta_all$EGFR_status <- ifelse(
  meta_all$patient_id %in% egfr_driver_patients,
  "EGFR_driver",
  "EGFR_WT"
)

# Keep primary tumor samples only, automatically LUAD/LUSC because of projects vector
meta_tumor <- meta_all %>% dplyr::filter(sample_type == "Primary Tumor")
counts_tumor <- counts_all[, meta_tumor$barcode]

# counts_tumor rownames are ENSG IDs, map them to HGNC symbols
res <- map_ensembl_to_hgnc(counts_tumor)

counts_tumor_mapped <- res$counts_mapped
gene_annotation     <- res$annotation

## 5) Normalize with DESeq2 (VST) ##############################################

dds <- DESeqDataSetFromMatrix(
  countData = counts_tumor_mapped,
  colData   = meta_tumor,
  design    = ~ EGFR_status
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

dds <- DESeq(dds)
vst_mat <- assay(vst(dds, blind = TRUE))

# Ensure EGFR present
if (!(egfr_gene %in% rownames(vst_mat))) {
  stop("EGFR disappeared after filtering – lower the filter threshold or check gene naming.")
}

egfr_expr <- vst_mat[egfr_gene, ]

## 6) Prepare gene lists (target + bypass + EGFR) ##############################

all_genes_of_interest <- unique(c(target_genes, bypass_genes, egfr_gene))
genes_present <- all_genes_of_interest[all_genes_of_interest %in% rownames(vst_mat)]
genes_missing <- setdiff(all_genes_of_interest, rownames(vst_mat))

if (length(genes_missing) > 0) {
  warning("These requested genes were not found in the expression matrix and will be skipped for relevant steps: ",
          paste(genes_missing, collapse = ", "))
}

target_genes_present <- target_genes[target_genes %in% rownames(vst_mat)]
bypass_genes_present <- bypass_genes[bypass_genes %in% rownames(vst_mat)]

## 7) Helper: safe correlation function #######################################

safe_cor <- function(x, y, method = "pearson") {
  if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
  if (var(x, na.rm = TRUE) == 0 || var(y, na.rm = TRUE) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = method, use = "pairwise.complete.obs"))
}

## 8) Per-gene analysis: co-expression + EGFR-driver association + LMs #########

meta_tumor$EGFR_status <- factor(meta_tumor$EGFR_status,
                                 levels = c("EGFR_WT","EGFR_driver"))
egfr_status     <- meta_tumor$EGFR_status
egfr_status_bin <- ifelse(egfr_status == "EGFR_driver", 1L, 0L)

summary_list <- list()

for (g in target_genes_present) {
  message("Processing target gene: ", g)
  g_expr <- vst_mat[g, ]
  
  # (A) Co-expression with EGFR
  pearson_egfr  <- safe_cor(g_expr, egfr_expr, method = "pearson")
  spearman_egfr <- safe_cor(g_expr, egfr_expr, method = "spearman")
  
  # (B) Co-expression with bypass genes (average correlation)
  bypass_corrs_pearson  <- c()
  bypass_corrs_spearman <- c()
  for (bg in bypass_genes_present) {
    bg_expr <- vst_mat[bg, ]
    bypass_corrs_pearson[bg]  <- safe_cor(g_expr, bg_expr, method = "pearson")
    bypass_corrs_spearman[bg] <- safe_cor(g_expr, bg_expr, method = "spearman")
  }
  
  mean_bypass_pearson  <- ifelse(length(bypass_corrs_pearson)  > 0,
                                 mean(bypass_corrs_pearson,  na.rm = TRUE),
                                 NA_real_)
  mean_bypass_spearman <- ifelse(length(bypass_corrs_spearman) > 0,
                                 mean(bypass_corrs_spearman, na.rm = TRUE),
                                 NA_real_)
  
  # (C) Association with EGFR-driver subtype (Wilcoxon + log2FC)
  expr_df <- data.frame(
    expr         = as.numeric(g_expr),
    EGFR_status  = egfr_status,
    EGFR_binary  = egfr_status_bin,
    project_id   = meta_tumor$project_id
  )
  
  wt_vals  <- expr_df$expr[expr_df$EGFR_status == "EGFR_WT"]
  drv_vals <- expr_df$expr[expr_df$EGFR_status == "EGFR_driver"]
  
  if (length(unique(na.omit(expr_df$EGFR_status))) == 2) {
    wtest <- wilcox.test(expr ~ EGFR_status, data = expr_df)
    wilcox_p <- wtest$p.value
  } else {
    wilcox_p <- NA_real_
  }
  
  mean_wt  <- mean(wt_vals,  na.rm = TRUE)
  mean_drv <- mean(drv_vals, na.rm = TRUE)
  log2fc_drv_vs_wt <- log2((mean_drv + 1e-6) / (mean_wt + 1e-6))
  
  ## (D) Linear model: expression ~ EGFR_status (+ project covariate) ##########
  lm_fit <- try(lm(expr ~ EGFR_status + project_id, data = expr_df), silent = TRUE)
  
  lm_beta_driver <- NA_real_
  lm_se_driver   <- NA_real_
  lm_p_driver    <- NA_real_
  lm_r2          <- NA_real_
  
  if (!inherits(lm_fit, "try-error")) {
    sm <- summary(lm_fit)
    if ("EGFR_statusEGFR_driver" %in% rownames(sm$coefficients)) {
      lm_beta_driver <- sm$coefficients["EGFR_statusEGFR_driver","Estimate"]
      lm_se_driver   <- sm$coefficients["EGFR_statusEGFR_driver","Std. Error"]
      lm_p_driver    <- sm$coefficients["EGFR_statusEGFR_driver","Pr(>|t|)"]
    }
    lm_r2 <- sm$r.squared
  }
  
  ## (E) Logistic model: EGFR_status_bin ~ expression (+ project_id) ###########
  glm_fit <- try(
    glm(EGFR_binary ~ expr + project_id, data = expr_df, family = binomial),
    silent = TRUE
  )
  
  glm_beta_expr <- NA_real_
  glm_se_expr   <- NA_real_
  glm_p_expr    <- NA_real_
  glm_or_expr   <- NA_real_
  
  if (!inherits(glm_fit, "try-error")) {
    smg <- summary(glm_fit)
    if ("expr" %in% rownames(smg$coefficients)) {
      glm_beta_expr <- smg$coefficients["expr","Estimate"]
      glm_se_expr   <- smg$coefficients["expr","Std. Error"]
      glm_p_expr    <- smg$coefficients["expr","Pr(>|z|)"]
      glm_or_expr   <- exp(glm_beta_expr)
    }
  }
  
  ## Store summary row for this gene ###########################################
  
  summary_list[[g]] <- tibble(
    gene                      = g,
    pearson_with_EGFR         = pearson_egfr,
    spearman_with_EGFR        = spearman_egfr,
    mean_bypass_pearson       = mean_bypass_pearson,
    mean_bypass_spearman      = mean_bypass_spearman,
    mean_expr_WT              = mean_wt,
    mean_expr_EGFR_driver     = mean_drv,
    log2FC_driver_vs_WT       = log2fc_drv_vs_wt,
    wilcox_p_driver_vs_WT     = wilcox_p,
    # Linear model metrics
    lm_beta_driver            = lm_beta_driver,
    lm_se_driver              = lm_se_driver,
    lm_p_driver               = lm_p_driver,
    lm_r2                     = lm_r2,
    # Logistic model metrics
    glm_beta_expr             = glm_beta_expr,
    glm_se_expr               = glm_se_expr,
    glm_p_expr                = glm_p_expr,
    glm_or_expr               = glm_or_expr
  )
  
  ## Save per-gene bypass correlation table ####################################
  
  if (length(bypass_corrs_pearson) > 0) {
    per_gene_bypass_df <- tibble(
      target_gene = g,
      bypass_gene = names(bypass_corrs_pearson),
      pearson     = as.numeric(bypass_corrs_pearson),
      spearman    = as.numeric(bypass_corrs_spearman)
    )
    write_csv(
      per_gene_bypass_df,
      file.path(output_dir, paste0("bypass_correlations_", g, ".csv"))
    )
  }
  
  ## 9) Plots for this gene ####################################################
  
  # (1) Scatter: EGFR vs target gene expression
  scatter_df <- data.frame(
    EGFR_expr   = egfr_expr,
    target_expr = g_expr,
    EGFR_status = egfr_status,
    project_id  = meta_tumor$project_id
  )
  
  p_scatter <- ggplot(scatter_df,
                      aes(x = EGFR_expr, y = target_expr, shape = EGFR_status)) +
    geom_point(alpha = 0.6) +
    theme_bw() +
    labs(
      title = paste0(g, " vs EGFR expression"),
      x     = "EGFR VST expression",
      y     = paste0(g, " VST expression")
    )
  
  ggsave(
    filename = file.path(plots_dir, paste0("scatter_EGFR_vs_", g, ".png")),
    plot     = p_scatter,
    width    = 5,
    height   = 4,
    dpi      = 300
  )
  
  # (2) Boxplot: target expression by EGFR-driver status
  p_box <- ggplot(expr_df, aes(x = EGFR_status, y = expr)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_bw() +
    labs(
      title     = paste0(g, " expression in EGFR-driver vs WT tumors"),
      x         = "",
      y         = paste0(g, " VST expression"),
      subtitle  = paste0("Wilcoxon p = ",
                         ifelse(is.na(wilcox_p), "NA", signif(wilcox_p, 3)))
    )
  
  ggsave(
    filename = file.path(plots_dir, paste0("boxplot_", g, "_EGFR_driver_vs_WT.png")),
    plot     = p_box,
    width    = 4,
    height   = 4,
    dpi      = 300
  )
}

## 10) Global summary table ####################################################

summary_df <- bind_rows(summary_list)

# Adjust p-values for multiple testing (BH)
summary_df <- summary_df %>%
  mutate(
    wilcox_p_adj_BH = p.adjust(wilcox_p_driver_vs_WT, method = "BH"),
    lm_p_adj_BH     = p.adjust(lm_p_driver,              method = "BH"),
    glm_p_adj_BH    = p.adjust(glm_p_expr,               method = "BH")
  )

write_csv(summary_df, file.path(output_dir, "target_gene_summary_with_LMs.csv"))

## 11) Global correlation heatmap (EGFR + targets + bypass) ####################

genes_for_heatmap <- unique(c(target_genes_present, bypass_genes_present, egfr_gene))
expr_mat_heat <- vst_mat[genes_for_heatmap, , drop = FALSE]

cor_mat <- cor(t(expr_mat_heat), use = "pairwise.complete.obs", method = "pearson")

cor_df <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("gene1","gene2","correlation")

p_heat <- ggplot(cor_df, aes(x = gene1, y = gene2, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1,1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title  = element_blank()
  ) +
  labs(title = "Pearson correlation (VST) among EGFR, targets, and bypass genes")

ggsave(
  filename = file.path(plots_dir, "correlation_heatmap_targets_bypass_EGFR.png"),
  plot     = p_heat,
  width    = 7,
  height   = 6,
  dpi      = 300
)

## 12) Volcano plot for EGFR-driver vs WT association ##########################
# Uses log2FC (EGFR_driver vs WT) on x-axis and -log10(Wilcoxon p) on y-axis
# Only the candidate target genes (target_genes_present) are included.

volcano_df <- summary_df %>%
  dplyr::filter(gene %in% target_genes_present) %>%
  mutate(
    neg_log10_p = -log10(wilcox_p_adj_BH),
    sig = ifelse(
      !is.na(wilcox_p_adj_BH) & (wilcox_p_adj_BH < 0.05) &
      !is.na(lm_p_adj_BH) & (lm_p_adj_BH < 0.1) &
      !is.na(glm_p_adj_BH) & (glm_p_adj_BH < 0.1),
      "significant",
      "not_significant"
    )
  )

# Handle infinite -log10(p) if p=0 (very rare)
if (any(!is.finite(volcano_df$neg_log10_p))) {
  max_finite <- max(volcano_df$neg_log10_p[is.finite(volcano_df$neg_log10_p)], na.rm = TRUE)
  volcano_df$neg_log10_p[!is.finite(volcano_df$neg_log10_p)] <- max_finite + 1
}

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FC_driver_vs_WT,
                        y = neg_log10_p,
                        label = gene,
                        shape = sig,
                        color = sig)) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "significant" = "red",
      "not_significant" = "black"
    )
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = subset(volcano_df, sig == "significant"),
    size = 2,
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    title = "Volcano plot: EGFR-driver vs WT (candidate targets)",
    x     = "log2(Expr_EGFR_driver / Expr_WT)",
    y     = "-log10(Wilcoxon fdr-adj-p-value)"
  )

ggsave(
  filename = file.path(plots_dir, "volcano_EGFR_driver_vs_WT_targets.png"),
  plot     = p_volcano,
  width    = 6,
  height   = 5,
  dpi      = 300
)

message("Done. Results written to: ", normalizePath(output_dir))
save.image(file = ".RData")
savehistory(file = ".Rhistory")
