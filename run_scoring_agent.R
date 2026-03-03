#!/usr/bin/env Rscript

# ============================================================
# Target Prioritization Agent (2-week MVP) - R
# Evidence sources (public): MyGene, Reactome, ClinicalTrials.gov, 
#         Europe PMC, ChEMBL
# Deterministic scoring via YAML rubrics + thresholds
# ============================================================

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(yaml)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(lubridate)
  library(readr)
  library(glue)
})

# -------------------------
# Utilities
# -------------------------

safe_req <- function(req, max_tries = 3, quiet = TRUE) {
  last_err <- NULL
  for (i in seq_len(max_tries)) {
    try({
      resp <- req |> req_perform()
      status <- resp_status(resp)
      if (status >= 200 && status < 300) return(resp)
      last_err <- paste("HTTP", status, resp_body_string(resp))
    }, silent = TRUE)
    Sys.sleep(0.5 * i)
  }
  if (!quiet) message("Request failed: ", last_err)
  return(NULL)
}

as_date_safe <- function(x) {
  # Attempts to parse common date strings
  if (is.null(x) || length(x) == 0 || all(is.na(x))) return(as.Date(NA))
  suppressWarnings({
    d <- ymd(x, quiet = TRUE)
    if (all(is.na(d))) d <- mdy(x, quiet = TRUE)
    if (all(is.na(d))) d <- as.Date(x)
  })
  d
}

clip <- function(x, lo, hi) pmax(lo, pmin(hi, x))

dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)

# -------------------------
# Config loaders
# -------------------------

load_yaml <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  yaml::read_yaml(path)
}

# -------------------------
# Target resolution (MyGene.info)
# -------------------------

resolve_targets_mygene <- function(target_symbols,
                                   override_path = "target_alias_overrides.csv",
                                   n_hits = 50) {
  base <- "https://mygene.info/v3/query"
  
  # optional overrides: columns input,canonical_symbol
  overrides <- tibble::tibble(input = character(), canonical_symbol = character())
  if (!is.null(override_path) && file.exists(override_path)) {
    overrides <- readr::read_csv(override_path, show_col_types = FALSE) |>
      dplyr::mutate(
        input = toupper(trimws(input)),
        canonical_symbol = toupper(trimws(canonical_symbol))
      ) |>
      dplyr::distinct(input, .keep_all = TRUE)
  }
  
  # helper to compute a deterministic "match score" for each hit
  score_hit <- function(input_u, hit) {
    sym <- toupper(hit$symbol %||% "")
    # alias may be absent, scalar, or vector
    alias_vec <- character(0)
    if (!is.null(hit$alias)) {
      alias_vec <- toupper(as.character(unlist(hit$alias)))
    }
    
    # ranking rules (bigger is better)
    if (sym == input_u) return(100L)                 # exact symbol match
    if (input_u %in% alias_vec) return(90L)          # exact alias match
    if (startsWith(sym, input_u)) return(40L)        # weak fallback
    if (any(startsWith(alias_vec, input_u))) return(30L)
    return(0L)
  }
  
  resolve_one <- function(sym_in) {
    input_u <- toupper(trimws(sym_in))
    
    # apply override if present
    q_term <- input_u
    if (nrow(overrides) > 0 && input_u %in% overrides$input) {
      q_term <- overrides$canonical_symbol[match(input_u, overrides$input)]
    }
    
    req <- httr2::request(base) |>
      httr2::req_url_query(
        q = q_term,
        species = "human",
        fields = "symbol,entrezgene,ensembl.gene,name,alias",
        size = n_hits
      )
    
    resp <- safe_req(req)
    if (is.null(resp)) {
      return(tibble::tibble(
        input = sym_in, symbol = NA_character_, entrez = NA_character_,
        ensembl = NA_character_, name = NA_character_, aliases = NA_character_
      ))
    }
    
    # parse without simplifyVector to avoid atomic-vector/data.frame pitfalls
    body <- jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
    hits <- body$hits %||% list()
    
    if (length(hits) == 0) {
      return(tibble::tibble(
        input = sym_in, symbol = NA_character_, entrez = NA_character_,
        ensembl = NA_character_, name = NA_character_, aliases = NA_character_
      ))
    }
    
    # score all candidates against ORIGINAL input (not overridden q_term)
    scores <- vapply(hits, function(h) score_hit(input_u, h), integer(1))
    best_i <- which.max(scores)
    best <- hits[[best_i]]
    
    # extract ensembl gene robustly
    ensg <- NA_character_
    if (!is.null(best$`ensembl.gene`)) {
      ensg <- as.character(unlist(best$`ensembl.gene`))[1]
    } else if (!is.null(best$ensembl) && !is.null(best$ensembl$gene)) {
      ensg <- as.character(unlist(best$ensembl$gene))[1]
    }
    
    aliases <- NA_character_
    if (!is.null(best$alias)) {
      aliases <- paste0(as.character(unlist(best$alias)), collapse = "|")
      if (aliases == "") aliases <- NA_character_
    }
    
    tibble::tibble(
      input = sym_in,
      symbol = as.character(best$symbol %||% NA_character_),
      entrez = as.character(best$entrezgene %||% NA_character_),
      ensembl = as.character(ensg),
      name = as.character(best$name %||% NA_character_),
      aliases = aliases
    )
  }
  
  purrr::map_dfr(target_symbols, resolve_one)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -------------------------
# Evidence retrieval
# -------------------------

# 1) Reactome pathways
get_reactome_pathways <- function(hgnc_symbol) {
  search_url <- "https://reactome.org/ContentService/search/query"
  
  req <- request(search_url) |>
    req_url_query(
      query = hgnc_symbol,
      species = "Homo sapiens",
      types = "ReferenceGeneProduct",
      cluster = "true"
    )
  
  resp <- safe_req(req)
  if (is.null(resp)) return(list(pathways = tibble(), n_pathways = 0))
  
  body <- resp |> resp_body_json(simplifyVector = TRUE)
  
  hits <- body$results %||% NULL
  if (is.null(hits) || length(hits) == 0) {
    return(list(pathways = tibble(stId = character(0), name = character(0)), n_pathways = 0))
  }
  
  # ---- Extract stIds robustly ----
  stIds <- character(0)
  
  if (is.data.frame(hits)) {
    if ("stId" %in% names(hits)) {
      stIds <- hits$stId
    } else if ("stId" %in% names(hits[[1]])) {
      # rare nested case
      stIds <- hits[[1]]$stId
    }
  } else if (is.list(hits)) {
    # list-of-lists
    stIds <- purrr::map_chr(hits, function(x) {
      if (is.list(x) && !is.null(x$stId)) return(as.character(x$stId))
      return(NA_character_)
    })
  }
  
  stIds <- unique(stIds)
  stIds <- stIds[!is.na(stIds) & stIds != ""]
  if (length(stIds) == 0) {
    return(list(pathways = tibble(stId = character(0), name = character(0)), n_pathways = 0))
  }
  
  # Keep only a few seed IDs to avoid lots of calls
  stIds <- stIds[1:min(3, length(stIds))]
  
  # ---- Fetch pathways for each stId ----
  pathways <- list()
  for (st in stIds) {
    pw_url <- glue("https://reactome.org/ContentService/data/pathways/low/diagram/entity/{st}")
    pw_req <- request(pw_url)
    pw_resp <- safe_req(pw_req)
    if (is.null(pw_resp)) next
    pw_body <- pw_resp |> resp_body_json(simplifyVector = TRUE)
    pathways <- c(pathways, list(pw_body))
  }
  
  # Flatten (pw_body is typically a list of pathway objects)
  flat <- list()
  for (item in pathways) {
    if (is.data.frame(item)) {
      flat <- c(flat, split(item, seq_len(nrow(item))))
    } else if (is.list(item)) {
      # item might already be list-of-pathways
      if (length(item) > 0 && is.list(item[[1]])) flat <- c(flat, item) else flat <- c(flat, list(item))
    }
  }
  
  if (length(flat) == 0) {
    df <- tibble(stId = character(0), name = character(0))
  } else {
    df <- purrr::map_dfr(flat, function(x) {
      if (is.data.frame(x)) x <- as.list(x[1, , drop = FALSE])
      if (!is.list(x)) x <- as.list(x)
      tibble(
        stId = as.character(x$stId %||% NA_character_),
        name = as.character(x$name %||% NA_character_)
      )
    }) |> distinct(stId, .keep_all = TRUE) |> filter(!is.na(stId))
  }
  
  list(pathways = df, n_pathways = nrow(df))
}

# 2) ClinicalTrials.gov
get_trials <- function(query, max_records = 50) {
  # v2 API: https://clinicaltrials.gov/api/v2/
  base <- "https://clinicaltrials.gov/api/v2/studies"
  
  req <- request(base) |>
    req_url_query(
      query.term = query,
      pageSize = max_records,
      fields = paste(c(
        "protocolSection.identificationModule.nctId",
        "protocolSection.identificationModule.briefTitle",
        "protocolSection.conditionsModule.conditions",
        "protocolSection.armsInterventionsModule.interventions",
        "protocolSection.statusModule.overallStatus",
        "protocolSection.statusModule.startDateStruct",
        "protocolSection.designModule.phases"
      ), collapse = ",")
    )
  
  resp <- safe_req(req)
  if (is.null(resp)) return(tibble())
  
  body <- resp |> resp_body_json(simplifyVector = TRUE)
  studies <- body$studies %||% list()
  if (length(studies) == 0) return(tibble())
  
  map_dfr(studies, function(st) {
    pm <- st$protocolSection %||% list()
    idm <- pm$identificationModule %||% list()
    cdm <- pm$conditionsModule %||% list()
    aim <- pm$armsInterventionsModule %||% list()
    sm <- pm$statusModule %||% list()
    dm <- pm$designModule %||% list()
    
    startDate <- sm$startDateStruct$date %||% NA_character_
    tibble(
      nct_id = idm$nctId %||% NA_character_,
      title = idm$briefTitle %||% NA_character_,
      conditions = paste0(cdm$conditions %||% character(0), collapse = "; "),
      interventions = paste0(map_chr(aim$interventions %||% list(), ~ .x$name %||% ""), collapse = "; "),
      status = sm$overallStatus %||% NA_character_,
      phase = paste0(dm$phases %||% character(0), collapse = "; "),
      start_date = as_date_safe(startDate)
    )
  }) |> distinct(nct_id, .keep_all = TRUE)
}

# 3) pubmed literature
get_pubmed <- function(query, retmax = 25, api_key = NULL) {
  base_search <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  base_summary <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
  
  # 1️⃣ Search PubMed
  req_search <- httr2::request(base_search) |>
    httr2::req_url_query(
      db = "pubmed",
      term = query,
      retmode = "json",
      retmax = retmax,
      api_key = api_key
    )
  
  resp_search <- safe_req(req_search)
  if (is.null(resp_search)) return(tibble())
  
  search_body <- jsonlite::fromJSON(httr2::resp_body_string(resp_search), simplifyVector = TRUE)
  
  ids <- search_body$esearchresult$idlist
  if (is.null(ids) || length(ids) == 0) return(tibble())
  
  # 2️⃣ Get metadata
  req_summary <- httr2::request(base_summary) |>
    httr2::req_url_query(
      db = "pubmed",
      id = paste(ids, collapse = ","),
      retmode = "json",
      api_key = api_key
    )
  
  resp_summary <- safe_req(req_summary)
  if (is.null(resp_summary)) return(tibble())
  
  summary_body <- jsonlite::fromJSON(httr2::resp_body_string(resp_summary), simplifyVector = FALSE)
  
  results <- summary_body$result
  uids <- results$uids
  
  if (is.null(uids) || length(uids) == 0) return(tibble())
  
  purrr::map_dfr(uids, function(uid) {
    r <- results[[uid]]
    
    pub_year <- NA_integer_
    if (!is.null(r$pubdate)) {
      pub_year <- suppressWarnings(as.integer(substr(r$pubdate, 1, 4)))
    }
    
    tibble(
      title = r$title %||% NA_character_,
      journal = r$fulljournalname %||% NA_character_,
      pub_year = pub_year,
      cited_by = NA_integer_,  # PubMed does not directly provide citation count
      pmid = uid,
      doi = if (!is.null(r$elocationid) && grepl("doi", r$elocationid, ignore.case = TRUE))
        r$elocationid else NA_character_
    )
  })
}

# 4) ChEMBL (public) - basic druggability proxy via target search + bioactivity count
get_chembl_target_summary <- function(hgnc_symbol) {
  # ChEMBL API
  # https://www.ebi.ac.uk/chembl/api/data/target/search?q=EGFR
  base <- "https://www.ebi.ac.uk/chembl/api/data/target/search"
  req <- request(base) |> req_url_query(q = hgnc_symbol, format = "json")
  resp <- safe_req(req)
  if (is.null(resp)) return(list(n_targets = 0, top_target = NULL))
  
  body <- resp |> resp_body_json(simplifyVector = TRUE)
  targets <- body$targets %||% list()
  if (length(targets) == 0) return(list(n_targets = 0, top_target = NULL))
  
  # pick first
  t0 <- targets[[1]]
  list(
    n_targets = length(targets),
    top_target = list(
      chembl_target_id = t0$target_chembl_id %||% NA_character_,
      pref_name = t0$pref_name %||% NA_character_,
      organism = t0$organism %||% NA_character_,
      target_type = t0$target_type %||% NA_character_
    )
  )
}

# -------------------------
# Evidence bundle builder
# -------------------------

build_evidence_bundle <- function(app, target_row) {
  sym <- target_row$symbol
  if (is.na(sym) || sym == "") sym <- target_row$input
  
  # Context queries (simple MVP)
  context_terms <- paste(app$context_terms, collapse = " OR ")
  disease_terms  <- paste(app$disease_terms, collapse = " OR ")
  therapy_terms  <- paste(app$therapy_terms, collapse = " OR ")
  
  # Trials: combine target + context
  trial_query <- glue('({sym}) AND ({disease_terms}) AND ({therapy_terms})')
  trials <- get_trials(trial_query, max_records = app$limits$trials_max)
  
  # Literature: target + context
  paper_query <- glue('{sym} AND ({disease_terms}) AND ({context_terms})')
  papers <- get_pubmed(
    query = paper_query,
    retmax = app$limits$papers_max
  )
  
  # Pathways
  reactome <- get_reactome_pathways(sym)
  
  # ChEMBL
  chembl <- get_chembl_target_summary(sym)
  
  # Derived metrics
  n_trials <- nrow(trials)
  n_trials_active <- sum(trials$status %in% c("RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"), na.rm = TRUE)
  most_recent_trial <- suppressWarnings(max(trials$start_date, na.rm = TRUE))
  if (is.infinite(most_recent_trial)) most_recent_trial <- as.Date(NA)
  
  n_papers <- nrow(papers)
  recent_papers <- papers |> filter(!is.na(pub_year) & pub_year >= (year(Sys.Date()) - 5))
  n_papers_5y <- nrow(recent_papers)
  
  bundle <- list(
    target = list(
      input = target_row$input,
      symbol = sym,
      entrez = target_row$entrez,
      ensembl = target_row$ensembl,
      name = target_row$name
    ),
    application = app,
    evidence = list(
      trials = trials,
      papers = papers,
      reactome_pathways = reactome$pathways,
      chembl = chembl
    ),
    metrics = list(
      n_trials = n_trials,
      n_trials_active = n_trials_active,
      most_recent_trial_start = as.character(most_recent_trial),
      n_papers = n_papers,
      n_papers_5y = n_papers_5y,
      n_pathways = reactome$n_pathways,
      chembl_targets_found = chembl$n_targets
    ),
    generated_at = as.character(Sys.time())
  )
  
  bundle
}

# -------------------------
# Scoring (deterministic, rubric-based)
# -------------------------

score_from_bins <- function(value, bins) {
  # bins: list of {min, score} sorted ascending by min
  # returns highest score whose min <= value
  s <- 0
  for (b in bins) {
    if (!is.null(b$min) && value >= b$min) s <- b$score
  }
  clip(s, 0, 5)
}

compute_scores <- function(bundle, scoring_cfg) {
  m <- bundle$metrics
  app <- bundle$application
  # convenience
  n_trials <- m$n_trials %||% 0
  n_trials_active <- m$n_trials_active %||% 0
  n_papers <- m$n_papers %||% 0
  n_papers_5y <- m$n_papers_5y %||% 0
  n_pathways <- m$n_pathways %||% 0
  chembl_found <- m$chembl_targets_found %||% 0
  
  # Category scoring using bins defined in scoring.yaml
  s <- list()
  
  # Biological rationale: pathways + literature in context
  s$biological_rationale <- round(mean(c(
    score_from_bins(n_pathways, scoring_cfg$biological_rationale$bins$n_pathways),
    score_from_bins(n_papers_5y, scoring_cfg$biological_rationale$bins$n_papers_5y)
  )), 1)
  
  # Genetic evidence: MVP placeholder (0 unless you plug TCGA/DepMap/RWD features)
  s$genetic_evidence <- scoring_cfg$genetic_evidence$default_score %||% 0
  
  # Therapeutic evidence: trials + context papers
  s$therapeutic_evidence <- round(mean(c(
    score_from_bins(n_trials, scoring_cfg$therapeutic_evidence$bins$n_trials),
    score_from_bins(n_trials_active, scoring_cfg$therapeutic_evidence$bins$n_trials_active)
  )), 1)
  
  # Druggability: ChEMBL target presence as proxy (improve later with target class)
  s$druggability <- score_from_bins(chembl_found, scoring_cfg$druggability$bins$chembl_targets_found)
  
  # Tooling: proxy via publications + trials activity + chembl presence
  s$tooling <- round(mean(c(
    score_from_bins(n_papers, scoring_cfg$tooling$bins$n_papers),
    score_from_bins(chembl_found, scoring_cfg$tooling$bins$chembl_targets_found)
  )), 1)
  
  # Translational feasibility: proxy via trials (esp active) + “broadness” of biology (pathways)
  s$translational_feasibility <- round(mean(c(
    score_from_bins(n_trials_active, scoring_cfg$translational_feasibility$bins$n_trials_active),
    score_from_bins(n_pathways, scoring_cfg$translational_feasibility$bins$n_pathways)
  )), 1)
  
  # Toxicity risk (higher=safer): MVP placeholder (3) until you add essentiality/GTEx/KO
  s$toxicity_risk <- scoring_cfg$toxicity_risk$default_score %||% 3
  
  # Competitive landscape (higher=better whitespace): invert trial crowding in this context
  # More trials => more crowded => lower score
  crowd_score <- score_from_bins(n_trials, scoring_cfg$competitive_landscape$bins$n_trials_crowding)
  # bins should be defined as "more trials => lower score", so we can use directly
  s$competitive_landscape <- crowd_score
  
  # Commercial viability: MVP proxy using (feasibility, whitespace, evidence strength)
  s$commercial_viability <- round(mean(c(
    s$translational_feasibility,
    s$competitive_landscape,
    s$therapeutic_evidence
  )), 1)
  
  # Weighted total
  weights <- scoring_cfg$weights
  # Ensure all keys exist
  cats <- names(weights)
  total <- 0
  for (cat in cats) {
    total <- total + (weights[[cat]] * (s[[cat]] %||% 0))
  }
  
  list(
    category_scores = s,
    total_weighted_score = round(total, 2),
    weights = weights
  )
}

# -------------------------
# Main runner
# -------------------------

run_agent <- function(app_path = "application.yaml",
                      scoring_path = "scoring.yaml",
                      targets_path = "targets.csv",
                      out_dir = "output") {
  
  dir_create(out_dir)
  dir_create(file.path(out_dir, "bundles"))
  
  app <- load_yaml(app_path)
  scoring_cfg <- load_yaml(scoring_path)
  
  targets <- read_csv(targets_path, show_col_types = FALSE) %>%
    mutate(target = as.character(target)) %>%
    filter(!is.na(target) & target != "") %>%
    distinct(target)
  
  resolved <- resolve_targets_mygene(targets$target)
  
  # Build bundles + score
  results <- pmap_dfr(resolved, function(input, symbol, entrez, ensembl, name, aliases = NA) {
    tr <- tibble(input = input, symbol = symbol, entrez = entrez, ensembl = ensembl, name = name, aliases = aliases)
    
    bundle <- build_evidence_bundle(app, tr)
    
    # Save evidence bundle JSON for auditability
    bundle_path <- file.path(out_dir, "bundles", paste0((bundle$target$symbol %||% bundle$target$input), ".json"))
    writeLines(toJSON(bundle, pretty = TRUE, auto_unbox = TRUE), bundle_path)
    
    score <- compute_scores(bundle, scoring_cfg)
    
    tibble(
      input = input,
      symbol = bundle$target$symbol %||% input,
      total_score = score$total_weighted_score,
      biological_rationale = score$category_scores$biological_rationale,
      genetic_evidence = score$category_scores$genetic_evidence,
      therapeutic_evidence = score$category_scores$therapeutic_evidence,
      druggability = score$category_scores$druggability,
      tooling = score$category_scores$tooling,
      translational_feasibility = score$category_scores$translational_feasibility,
      toxicity_risk = score$category_scores$toxicity_risk,
      competitive_landscape = score$category_scores$competitive_landscape,
      commercial_viability = score$category_scores$commercial_viability,
      evidence_bundle = bundle_path
    )
  })
  
  ranked <- results %>% arrange(desc(total_score), desc(therapeutic_evidence), desc(biological_rationale))
  out_csv <- file.path(out_dir, "ranked_targets.csv")
  write_csv(ranked, out_csv)
  
  # Also save a JSON summary
  out_json <- file.path(out_dir, "ranked_targets.json")
  writeLines(toJSON(ranked, pretty = TRUE, auto_unbox = TRUE), out_json)
  
  message("Done. Outputs:\n- ", out_csv, "\n- ", out_json, "\n- bundles in ", file.path(out_dir, "bundles"))
  invisible(ranked)
}

# -------------------------
# CLI entry (minimal)
# -------------------------
args <- commandArgs(trailingOnly = TRUE)
# Allow: Rscript run_agent.R --app application.yaml --scoring scoring.yaml --targets targets.csv --out output
arg_get <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) return(default)
  args[[idx + 1]]
}

app_path <- arg_get("--app", "application.yaml")
scoring_path <- arg_get("--scoring", "scoring.yaml")
targets_path <- arg_get("--targets", "targets.csv")
out_dir <- arg_get("--out", "output")

run_agent(app_path, scoring_path, targets_path, out_dir)