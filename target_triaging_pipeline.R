#################################################################################################################
# IN SILICO TARGET TRIAGING PIPELINE
# For EGFR-mutant NSCLC sensitization targets
#
# Pipeline Stages (NR is Non-Responder, R is Responder)
# 1. Input — ingest the 20+ genes and their methylation comparison results between NR and R groups 
# 2. Annotation & biological context — annotate each gene for function, pathway, expression in LUAD/NSCLC, 
#    EGFR interaction, subcellular localization, known resistance/response literature.
#    Drug / target annotation: ChEMBL, DrugBank, Open Targets, UniProt, PDB for structures.
#    Pathway / network: Reactome, KEGG, STRING for PPI networks and pathway embedding.
# 3. Evidence mining (multi-layer) — gather orthogonal evidence:
#     - clinical RWE signal strength (co-occurrence with progression or response)
#     - tumor genomics: TCGA-LUAD, AACR GENIE (for prevalence), cBioPortal for 
#               co-occurrence, COSMIC for mutation context
#     - functional genomics: DepMap (CRISPR knockouts, RNAi), Project SCORE, Achilles datasets, DRIVE
#     - Expression & proteomics: TCGA RNA-seq, CPTAC proteomics (if available), GTEx for normal tissue 
#               expression (toxicity risk).
#     - epigenetic evidence: methylation correlation with expression and outcome
#     - phospho / proteomic evidence if available
# 4. Druggability & tooling — structural druggability, ligandability, antibodies/assays available, tool compounds, 
#    tractability class (small molecule vs biologic). 
#    Data sources: commercial suppliers (CST, Abcam, Sigma, Selleck), BenchSci for reagent evidence.
# 5. Competitive landscape & translational feasibility — existing drugs/clinical trials, patents, prevalence 
#    in patient population, biomarker assay readiness.
#    Data sources: ClinicalTrials.gov, patent databases (Lens.org, Google Patents)
# 6. Risk & technical feasibility — gene essentiality (toxicity risk), redundancy/compensation in pathways, 
#    availability of selective modulators, assay difficulty.
# 7. Scoring & rank ordering — compute normalized subscores, aggregate (weighted), produce rank and uncertainty metrics.
# 8. Sensitivity analysis & robustness — test different weighting schemes, bootstrap subscores to measure rank stability.
# 9. Output & recommendation — prioritized list with rationale & recommended next experiments for top N (e.g., top 3–5).
#################################################################################################################

# Load packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(purrr)
  library(ggplot2)
  library(scales)
})

###############################################
# 1. Prepare Input Data
###############################################
#
# 1. generate drug_target_annotations.csv which is the “one-stop” table for target tractability and tools
# it should contain:
#     targets' core IDs, target class & biology, druggability/tractability, 
###############################################
# 2. Load Input Data
###############################################

# User provides these CSV/TXT files
guardant_file      <- "pathway_genes_baseline_methyl_nofilt_NRvR_2025Q3.txt"      # potential targets
depmap_crispr_file <- "/domino/datasets/local/cohort_comparisons/DepMap/25Q3/CRISPRGeneEffect.csv"      # DepMap CERES/Chronos
depmap_rnai_file   <- "domino/datasets/local/cohort_comparisons/DepMap/D2_combined_gene_dep_scores.csv" # DepMap DEMETER2
depmap_mutants_file<- "depmap_EGFR_mutant_lines.csv"  # List of EGFR-mutant cell lines
tcga_prevalence    <- "tcga_luad_prevalence.csv"      # query TCGA NSCLC samples to get the frequency of promoter or 
                                                      # enhancer methylation events for each gene
drug_target_file   <- "drug_target_annotations.csv"
assay_availability <- "assay_annotations.csv"
mouse_knockout     <- "mouse_KO_phenotypes.csv"
normal_expression  <- "gtex_normal_expression.csv"

# Load
guardant   <- read_csv(guardant_file)
depmap     <- read_csv(depmap_file)
depmap_mut <- read_csv(depmap_mutants_file)
tcga       <- read_csv(tcga_prevalence)
drugdb     <- read_csv(drug_target_file)
assaydb    <- read_csv(assay_availability)
mouseKO    <- read_csv(mouse_knockout)
gtex       <- read_csv(normal_expression)

###############################################
# 2. Helper Normalization Functions
###############################################

# Percentile normalization 0–1
norm01 <- function(x) {
  if (length(unique(x)) == 1) return(rep(0.5, length(x)))
  return((rank(x, na.last = "keep") - 1) / (sum(!is.na(x)) - 1))
}

# Inverse percentile (higher = worse)
inv_norm01 <- function(x){
  out <- norm01(x)
  return(1 - out)
}

###############################################
# 3. Compute Subscores
###############################################

df <- guardant %>%
  rename(gene = Gene)

# 3.1 CRWE subscore
df <- df %>%
  mutate(
    s_CRWE = 0.5 * norm01(pct_nonresponders) +
      0.3 * norm01(log2_vaf_ratio) +
      0.2 * ifelse(methylation_significant == 1, 1, 0)
  )

# 3.2 Functional genomics (DepMap)
depmap_EGFR <- depmap %>% 
  filter(cell_line %in% depmap_mut$cell_line)

fun_scores <- depmap_EGFR %>%
  group_by(gene) %>%
  summarize(dep_mean = mean(gene_effect, na.rm = TRUE)) %>%
  mutate(s_FUN = inv_norm01(dep_mean))   # more negative gene effect = more essential

df <- df %>%
  left_join(fun_scores, by = "gene")

# 3.3 Biological relevance score placeholder (user may add pathway data)
# For now use curated annotations from your drug_target_file
df <- df %>%
  left_join(select(drugdb, gene, pathway_score, egfr_distance_score), by = "gene") %>%
  mutate(
    s_BIO = 0.6 * norm01(pathway_score) +
      0.4 * inv_norm01(egfr_distance_score)   # closer to EGFR = stronger
  )

# 3.4 TCGA prevalence
df <- df %>%
  left_join(tcga, by = "gene") %>%
  mutate(s_PREV = norm01(prevalence))

# 3.5 Druggability
df <- df %>%
  left_join(drugdb, by = "gene") %>%
  mutate(
    s_DRUG = 0.5 * norm01(druggability_score) +
      0.5 * norm01(tool_compound_score)
  )

# 3.6 Assay readiness
df <- df %>%
  left_join(assaydb, by = "gene") %>%
  mutate(
    s_ASSAY = 0.5 * norm01(antibody_quality) +
      0.5 * norm01(assay_complexity_score)
  )

# 3.7 Competitive landscape (optional sign flip)
df <- df %>%
  mutate(s_COMP = inv_norm01(competition_intensity))

###############################################
# 4. Toxicity Risk Subscore
###############################################

# 4.1 Normal tissue expression (higher = higher toxicity risk)
tox_expr <- gtex %>%
  group_by(gene) %>%
  summarize(normal_expr = mean(tpm, na.rm = TRUE)) %>%
  mutate(expr_risk = norm01(normal_expr))

# 4.2 Mouse knockout lethality
mouse_risk <- mouseKO %>%
  mutate(KO_risk = case_when(
    phenotype == "embryonic lethal" ~ 1,
    phenotype == "postnatal lethal" ~ 0.8,
    phenotype == "subfertile" ~ 0.4,
    TRUE ~ 0.1
  ))

# 4.3 DepMap essentiality as toxicity component
dep_essential <- depmap %>%
  group_by(gene) %>%
  summarize(pan_dep = mean(gene_effect, na.rm = TRUE)) %>%
  mutate(essential_risk = inv_norm01(pan_dep))

# Combine toxicity components
tox <- tox_expr %>%
  left_join(mouse_risk, by = "gene") %>%
  left_join(dep_essential, by = "gene") %>%
  mutate(
    KO_risk = ifelse(is.na(KO_risk), 0.2, KO_risk),
    toxicity_risk = 0.4 * expr_risk +
      0.3 * KO_risk +
      0.3 * essential_risk
  ) %>%
  select(gene, toxicity_risk)

df <- df %>% left_join(tox, by = "gene")

###############################################
# 5. Composite Score
###############################################

weights <- list(
  CRWE = 0.20,
  BIO  = 0.18,
  FUN  = 0.16,
  DRUG = 0.14,
  PREV = 0.10,
  ASSAY= 0.08,
  COMP = 0.06,
  RISK = 0.08   # penalty weight
)

df <- df %>%
  mutate(
    composite_score =
      weights$CRWE * s_CRWE +
      weights$BIO  * s_BIO  +
      weights$FUN  * s_FUN  +
      weights$DRUG * s_DRUG +
      weights$PREV * s_PREV +
      weights$ASSAY* s_ASSAY +
      weights$COMP * s_COMP -
      weights$RISK * toxicity_risk
  ) %>%
  arrange(desc(composite_score)) %>%
  mutate(rank = row_number())

###############################################
# 6. Output
###############################################

write_csv(df, "gene_ranking_output.csv")

print("Ranking complete. Output saved to gene_ranking_output.csv.")

###############################################
# 7. Optional: Rank Stability via Bootstrap
###############################################

bootstrap_rank <- function(df, B = 500){
  genes <- df$gene
  ranks_matrix <- matrix(0, nrow = length(genes), ncol = B)
  rownames(ranks_matrix) <- genes
  
  for (b in 1:B) {
    # randomly perturb weights ±20%
    perturbed <- weights
    for (w in names(perturbed)) {
      perturbed[[w]] <- perturbed[[w]] * runif(1, 0.8, 1.2)
    }
    # recompute composite
    score <- df %>%
      mutate(score =
               perturbed$CRWE * s_CRWE +
               perturbed$BIO  * s_BIO  +
               perturbed$FUN  * s_FUN  +
               perturbed$DRUG * s_DRUG +
               perturbed$PREV * s_PREV +
               perturbed$ASSAY* s_ASSAY +
               perturbed$COMP * s_COMP -
               perturbed$RISK * toxicity_risk) %>%
      arrange(desc(score))
    
    ranks_matrix[score$gene, b] <- 1:nrow(score)
  }
  return(ranks_matrix)
}

# To run:
# ranks_mat <- bootstrap_rank(df, B = 300)
# saveRDS(ranks_mat, "bootstrap_ranks.rds")

###############################################
# END OF PIPELINE
###############################################
