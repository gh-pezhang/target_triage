###########################################################################
## script to evaluate how dependent EGFR-mutant NSCLC cell lines are on
## each of the genes in the target gene list
## workflow that uses DepMap data to:
## 1. Load CRISPR (or RNAi) dependency scores, cell line metadata, and mutation data
## 2. Define EGFR-driver NSCLC cell lines (lung + EGFR oncogenic mutation)
## 3. Extract dependency scores for target genes
## 4. Compare EGFR-mutant vs EGFR-WT lung lines
## 5. Generate summary tables and plots (boxplots + heatmap)
###########################################################################

# 0. User parameters ###########################################
setwd("/mnt/target_triage")

# Output directories
output_dir <- "./osi_resistance_targets_depmap_results"
plots_dir  <- file.path(output_dir, "plots")
if (!dir.exists(output_dir)) dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(plots_dir)) dir.create(plots_dir,  showWarnings = FALSE, recursive = TRUE)

target_genes <- c(
  "IL16", "PTPN6", "IL1A", "PLAT", "LY96", "PTH1R", "ID3", "CD7", 
  "PTEN", "CD47", "GALM", "IL18R1", "RHOH", "CARTPT", "IL2RB", 
  "ADORA2B", "NT5E", "ICOSLG", "MET", "CDKN1C", "TNFRSF12A", 
  "PFKP", "AXL", "HK2", "CELSR2", "DYRK3", "CYP39A1", "LTBR", "HDAC9"
)

# 1. Load packages and DepMap tables
if (!requireNamespace("ggforce", quietly = TRUE)) install.packages("ggforce")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(pheatmap)
  library(ggrepel)
  library(gridExtra)
  library(ggforce)
})

######### helper function #########
# Paginate a facet_wrap plot (3x5 per page) and return all pages
# X is a dataframe with samples in rows and features in columns
plot_feature_boxplots <- function(
    X, y, features = NULL, ylab, rawp = "N", title = NULL, outfile = NULL,
    rows = 3, cols = 5, logy = "N"
) {
  stopifnot(nrow(X) == length(y))
  y <- factor(y)
  if (nlevels(y) != 2) stop("y must have exactly two classes.")
  classes <- levels(y)
  
  # choose 10 features if not provided
  if (is.null(features)) features <- colnames(X)[seq_len(min(10, ncol(X)))]
  
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(ggplot2); library(stringr); library(ggforce)
  })
  
  # tidy data
  df_long <- as.data.frame(X[, features, drop = FALSE]) %>%
    mutate(.y = y) %>%
    pivot_longer(cols = all_of(features), names_to = "feature", values_to = "value") %>%
    mutate(feature = factor(feature, levels = features),
           .y = factor(.y, levels = classes))
  
  # Wilcoxon p-values per feature, then FDR adjust
  pvals <- df_long %>%
    group_by(feature) %>%
    summarise(
      p = {
        v1 <- value[.y == classes[1]]
        v2 <- value[.y == classes[2]]
        tryCatch(wilcox.test(v1, v2, exact = FALSE)$p.value, error = function(e) NA_real_)
      },
      .groups = "drop"
    ) %>%
    mutate(
      p_fdr = p.adjust(p, method = "fdr", n = ncol(X)),
      p_lab = dplyr::case_when(
        rawp == "N" ~ paste0("FDR p = ", formatC(p_fdr, format = "e", digits = 2)),
        TRUE        ~ paste0("RAW p = ", formatC(p,     format = "e", digits = 2))
      )
    )
  
  # y-position for annotation
  y_tops <- df_long %>%
    group_by(feature) %>%
    summarise(ymax = max(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(y_pos = ymax + 0.05 * (abs(ymax) + 1e-9))
  
  ann <- left_join(pvals, y_tops, by = "feature")
  
  # palette for two classes
  pal <- c("cyan", "red")
  
  # Base (unfaceted) plot
  # Add a small offset so zeros can be plotted
  min_pos <- min(df_long$value[df_long$value > 0], na.rm = TRUE)
  eps <- if (is.finite(min_pos)) min_pos / 10 else 1e-8   # fallback if all zeros
  
  # Use the offset in BOTH the data and your annotations:
  df_long <- df_long %>%
    mutate(value_plot = value + eps)
  
  ann <- ann %>%
    mutate(y_pos_plot = y_pos + eps)  # keep annotation above boxes on the log scale
  
  base_plot <- ggplot(df_long, aes(x = .y, y = value_plot, fill = .y)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.6) +
    scale_fill_manual(values = pal, name = "Class") +
    geom_text(
      data = ann, aes(x = 1.5, y = y_pos_plot, label = p_lab),
      inherit.aes = FALSE, size = 3
    ) +
    labs(title = title, x = NULL, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(size = 8, face = "bold", hjust = 0.5)
    )
  
  # add log-scale option
  if (toupper(logy) == "Y") {
    base_plot <- base_plot + scale_y_log10()
  }
  
  # Paginate
  paged_proto <- base_plot +
    ggforce::facet_wrap_paginate(~ feature, nrow = rows, ncol = cols, scales = "free_y")
  
  total_pages <- ggforce::n_pages(paged_proto)
  
  pages <- vector("list", total_pages)
  for (i in seq_len(total_pages)) {
    pg <- base_plot +
      ggforce::facet_wrap_paginate(
        ~ feature, nrow = rows, ncol = cols, page = i, scales = "free_y"
      ) +
      labs(subtitle = sprintf("Page %d / %d", i, total_pages))
    pages[[i]] <- pg
  }
  
  # Optional saving
  if (!is.null(outfile)) {
    if (grepl("\\.pdf$", outfile, ignore.case = TRUE)) {
      pdf(outfile, width = 12, height = 9, onefile = TRUE)
      on.exit(dev.off(), add = TRUE)
      for (pg in pages) print(pg)
    } else {
      parts <- tools::file_path_sans_ext(outfile)
      ext   <- paste0(".", tools::file_ext(outfile))
      if (ext == ".") ext <- ".png"
      for (i in seq_along(pages)) {
        f_i <- sprintf("%s_p%02d%s", parts, i, ext)
        ggsave(filename = f_i, plot = pages[[i]], width = 12, height = 9, units = "in", dpi = 300)
      }
    }
  }
  
  return(pages)
}

# Paths to DepMap files (adjust to your filenames)
yearqrt <- '25Q3'
depmap_dir      <- paste("/domino/datasets/local/cohort_comparisons/DepMap/", yearqrt, sep="")

file_dep       <- file.path(depmap_dir, "CRISPRGeneEffect.csv")       # gene dependency
file_meta      <- file.path(depmap_dir, "Model.csv")                  # cell line metadata
file_mut       <- file.path(depmap_dir, "OmicsSomaticMutations.csv")  # mutation calls

# Read in tables
dep_df  <- read_csv(file_dep)
colnames(dep_df)[1] <- 'ModelID'
meta_df <- read_csv(file_meta)
mut_df  <- read_csv(file_mut)

# 2. Restrict to NSCLC (lung) cell lines
# Inspect meta to see how lung is encoded
table(meta_df$OncotreeLineage)
table(meta_df$OncotreePrimaryDisease)

# Define NSCLC / lung subset
nsclc_meta <- meta_df %>%
  filter(
    grepl("lung", OncotreeLineage, ignore.case = TRUE) &
      grepl("Non-Small Cell Lung Cancer", OncotreePrimaryDisease, ignore.case = TRUE)
  )

nsclc_ids <- nsclc_meta$ModelID
length(nsclc_ids)

# 3. Define EGFR-driver NSCLC lines from mutation data
# limit to nsclc line mutations & dependencies
mut_df <- mut_df %>% filter(ModelID %in% nsclc_ids)
dep_df <- dep_df %>% filter(ModelID %in% nsclc_ids)

# Look at mutation columns
colnames(mut_df)

# Try to identify "impact" / "oncogenic" annotation columns
# Example relevant columns might be: 'Variant_annotation', 'isTCGAhotspot', 'OncoKB_Annotation', etc.

egfr_mut_dmap <- mut_df %>%
  filter(HugoSymbol == "EGFR")

# If there is an "OncoKB_Annotation" or similar oncogenic column:
#oncogenic_col <- intersect(
#  c("OncoKB_Annotation","onco_annotation","ONCOGENIC"), 
#  colnames(egfr_mut)
#)
oc_col <- 'VepClinSig'

if (!is.na(oc_col)) {
  egfr_driver_dmap <- egfr_mut_dmap %>%
    filter( grepl("pathogenic", .data[[oc_col]]))
} else {
  # Fallback: use simple variant patterns for activating EGFR mutations
  egfr_driver_dmap <- egfr_mut_dmap %>%
    filter(
      VariantType == "SNV" &
        (
          grepl("L858R", ProteinChange) |
            grepl("T790M", ProteinChange) |
            grepl("G719[ACS]", ProteinChange) |
            grepl("L861Q", ProteinChange) |
            grepl("S768I", ProteinChange) |
            grepl("E709", ProteinChange)
        ) |
      grepl("deletion", VariantType) & grepl("inframe_deletion", MolecularConsequence) |
      grepl("insertion", VariantType)
    )
}

egfr_driver_nsclc_ids <- unique(egfr_driver_dmap$ModelID)
length(egfr_driver_nsclc_ids)

# EGFR-WT NSCLC lines
egfr_wt_nsclc_ids <- setdiff(nsclc_ids, egfr_driver_nsclc_ids)
length(egfr_wt_nsclc_ids)

# 4. Reshape dependency data and focus on lung NSCLC
# Ensure ModelID column is present
colnames(dep_df)[1:5]

# Keep only target genes + ModelID
colnames(dep_df) <- gsub(" \\(\\d+\\)$", "", colnames(dep_df))  #take out the GeneID in colnames
genes_for_dep <- intersect(target_genes, colnames(dep_df))
length(genes_for_dep)

dep_df_long <- dep_df %>%
  select(ModelID, all_of(genes_for_dep)) %>%
  tidyr::pivot_longer(
    cols = -ModelID,
    names_to = "gene",
    values_to = "dependency"
  )

# Add EGFR-driver status and lineage info
nsclc_meta_sub <- nsclc_meta %>%
  select(ModelID, cell_line_name = CellLineName, lineage = OncotreeLineage, primary_disease = OncotreePrimaryDisease)

dep_df_long <- dep_df_long %>%
  left_join(nsclc_meta_sub, by = "ModelID") %>%
  mutate(
    EGFR_status = ifelse(ModelID %in% egfr_driver_nsclc_ids,
                         "EGFR_driver",
                         "EGFR_WT")
  )

# 5. Summary statistics + wilcoxon test per gene
summary_dep <- dep_df_long %>%
  group_by(gene, EGFR_status) %>%
  summarize(
    n        = n(),
    median_dep = median(dependency, na.rm = TRUE),
    sd_dep   = sd(dependency,   na.rm = TRUE),
    .groups = "drop"
  )

# Wide summary: EGFR_driver vs EGFR_WT
summary_dep_wide <- summary_dep %>%
  tidyr::pivot_wider(
    id_cols = gene,
    names_from = EGFR_status,
    values_from = c(n, median_dep, sd_dep),
    names_sep = "_"
  )

# Per-gene wilconxon test
wilcoxon_list <- dep_df_long %>%
  group_by(gene) %>%
  summarize(
    wilcoxon_p = tryCatch(
      wilcox.test(dependency ~ EGFR_status, alternative = "two.sided")$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

summary_dep_final <- summary_dep_wide %>%
  left_join(wilcoxon_list, by = "gene") %>%
  mutate(
    # Direction: more negative = stronger dependency
    perc_delta_median_dep = (median_dep_EGFR_driver - median_dep_EGFR_WT)/abs(median_dep_EGFR_WT),
    wilcoxon_p_adj_BH = p.adjust(wilcoxon_p, method = "BH")
  )

# Save summary
write_csv(summary_dep_final,
          file.path(output_dir,
                    "target_dependency_summary.csv"))

# 6. Boxplots per gene (EGFR-driver vs WT)
dep_df_wide <- dep_df_long %>%
  select(ModelID, EGFR_status, gene, dependency) %>%
  tidyr::pivot_wider(
    names_from  = gene,
    values_from = dependency
  )

pdf_file <- paste(plots_dir, "/target_gene_dependency_boxplots_depmap_", yearqrt, ".pdf", sep="")
p <- plot_feature_boxplots(dep_df_wide%>%select(-c("ModelID", "EGFR_status")), dep_df_wide$EGFR_status, 
                           features = genes_for_dep, ylab="CRISPR dependency (lower = more essential)", 
                           title="Gene Dependency in NSCLC cell lines", rawp = "Y",
                           outfile=pdf_file,
                           logy="N")

message("Multi-page PDF written to: ", pdf_file)


# 7. Heatmap of median dependency (genes × EGFR-status)
median_dep_by_status <- dep_df_long %>%
  group_by(gene, EGFR_status) %>%
  summarize(median_dep = median(dependency, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = EGFR_status,
    values_from = median_dep
  )

median_dep_mat <- median_dep_by_status %>%
  column_to_rownames("gene") %>%
  as.matrix()

h <- pheatmap(median_dep_mat,
         cluster_rows = TRUE,
         main = "median dependency by EGFR status (NSCLC)")
ggsave(file.path(plots_dir, "heatmap_median_dep_by_EGFR_status.png"),
       plot = h, width = 4, height = 5, dpi = 300)

# Volcano-style plot for dependency differences
volcano_dep <- summary_dep_final %>%
  mutate(
    neg_log10_p = -log10(wilcoxon_p),
    delta_dep = ifelse(
      #!is.na(wilcoxon_p) &
      #  wilcoxon_p < 0.1 &      # softer threshold
        perc_delta_median_dep < -0.25,
      ">25% more dependent in EGFR_driver",
      "not >25% more dependent in EGFR_driver"
    )
  )

v <- ggplot(volcano_dep,
       aes(x = perc_delta_median_dep, y = neg_log10_p, color = delta_dep, label = gene)) +
  geom_point() +
  scale_color_manual(values = c(">25% more dependent in EGFR_driver" = "red", 
                                "not >25% more dependent in EGFR_driver" = "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_text_repel(
    data = subset(volcano_dep, delta_dep == ">25% more dependent in EGFR_driver"),
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    title = "EGFR-driver vs WT NSCLC: target gene dependency",
    x = "median dependency (EGFR_driver − WT)/abs(WT)",
    y = "-log10 p-value"
  )

ggsave(file.path(plots_dir, "volcano_plot_dependency_EGFR_driver_vs_WT.png"),
       plot = v, width = 10, height = 5, dpi = 300)

save.image(file="depmap.RData")
