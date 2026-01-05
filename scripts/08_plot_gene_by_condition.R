## ============================================================
## 08_plot_gene_by_condition.R
## Gene-level visualization per condition, using Script 03 objects:
##   - raw mean per time (replicates averaged within time×GEAR)
##   - per-dataset facets
##   - OPTIONAL: min–max normalization (0..1) per GEAR with gray overlays
##              and mean across datasets in black
##
## NEW (optional):
##   - attach promoter motif counts for the gene (if motif tables exist)
##   - attach a defensible rhythmicity "strength" summary:
##       min/median adjusted p-values (padj) per method across SUPPORTING gears
##
## REQUIREMENTS:
##   Source Script 03 first (expr_* and annot_* objects exist)
##
## USAGE:
##   res <- plot_gene_by_condition("AT5G61380", "12L12D", mode="raw")
##   res$plot_all; res$plot_panels
##   res$motif_summary
##   res$rhythmicity_strength
##
##   resn <- plot_gene_by_condition("AT5G61380", "12L12D", mode="norm01")
##   resn$plot_norm_overlay
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## ----------------------------
## 1) Condition object getter
## ----------------------------

get_condition_data <- function(condition) {
  cond <- match.arg(condition, choices = c("12L12D", "12L12D_LL", "16L8D", "8L16D"))
  
  if (cond == "12L12D") {
    if (!exists("expr_12L12D") || !exists("annot_12L12D")) {
      stop("expr_12L12D / annot_12L12D not found. Run 03_load_PPRclock_all_conditions.R first.")
    }
    return(list(expr = expr_12L12D, annot = annot_12L12D, label = "LD (12L12D)"))
  }
  
  if (cond == "12L12D_LL") {
    if (!exists("expr_12L12D_LL") || !exists("annot_12L12D_LL")) {
      stop("expr_12L12D_LL / annot_12L12D_LL not found. Run 03_load_PPRclock_all_conditions.R first.")
    }
    return(list(expr = expr_12L12D_LL, annot = annot_12L12D_LL, label = "LD→LL (12L12D→LL)"))
  }
  
  if (cond == "16L8D") {
    if (!exists("expr_16L8D") || !exists("annot_16L8D")) {
      stop("expr_16L8D / annot_16L8D not found. Run 03_load_PPRclock_all_conditions.R first.")
    }
    return(list(expr = expr_16L8D, annot = annot_16L8D, label = "Long Day (16L8D)"))
  }
  
  if (cond == "8L16D") {
    if (!exists("expr_8L16D") || !exists("annot_8L16D")) {
      stop("expr_8L16D / annot_8L16D not found. Run 03_load_PPRclock_all_conditions.R first.")
    }
    return(list(expr = expr_8L16D, annot = annot_8L16D, label = "Short Day (8L16D)"))
  }
}

## ----------------------------
## 2) Input checks + data builder
## ----------------------------

.validate_inputs <- function(gene_id, expr_mat, annot) {
  if (!is.matrix(expr_mat) && !is.data.frame(expr_mat)) stop("expr_mat must be a matrix or data.frame.")
  if (!"column_name" %in% colnames(annot)) stop("annot must contain 'column_name'.")
  if (!all(c("gear_id", "time") %in% colnames(annot))) stop("annot must contain 'gear_id' and 'time'.")
  if (!gene_id %in% rownames(expr_mat)) stop("Gene not found in expr_mat: ", gene_id)
  if (!all(annot$column_name %in% colnames(expr_mat))) {
    bad <- setdiff(annot$column_name, colnames(expr_mat))
    stop("These annot$column_name are missing in expr_mat columns: ", paste(head(bad, 10), collapse = ", "),
         if (length(bad) > 10) " ..." else "")
  }
  invisible(TRUE)
}

.clean_agi <- function(x) toupper(trimws(sub("\\..*$", "", x)))

.build_gene_long <- function(gene_id, expr_mat, annot) {
  gene_values <- expr_mat[gene_id, ]
  data.frame(
    AGI         = gene_id,
    expression  = as.numeric(gene_values),
    column_name = names(gene_values),
    stringsAsFactors = FALSE
  ) %>%
    left_join(annot, by = "column_name") %>%
    filter(is.finite(time), is.finite(gear_id))
}

## Replicates averaged within each time×GEAR
.summarise_by_gear_time <- function(df_long) {
  df_long %>%
    group_by(gear_id, time) %>%
    summarise(
      expression_mean = mean(expression, na.rm = TRUE),
      n_reps = n(),
      .groups = "drop"
    ) %>%
    arrange(gear_id, time)
}

## ----------------------------
## 2.5) NEW: Motifs + rhythmicity strength helpers
## ----------------------------

.safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  # Use data.table fread if available (faster / safer for big tables)
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::fread(path, data.table = FALSE))
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.get_motif_summary <- function(gene_id,
                               condition,
                               motif_dir = "tables/motif_scan_per_photoperiod") {
  gene_id <- .clean_agi(gene_id)
  motif_file <- file.path(motif_dir, paste0("motif_counts_", condition, ".csv"))
  mot <- .safe_read_csv(motif_file)
  if (is.null(mot) || !"AGI" %in% colnames(mot)) return(NULL)
  
  mot$AGI <- .clean_agi(mot$AGI)
  out <- mot[mot$AGI == gene_id, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  
  # Collapse duplicates defensively (should already be 1 row per gene)
  motif_cols <- setdiff(colnames(out), "AGI")
  out2 <- data.frame(AGI = gene_id, stringsAsFactors = FALSE)
  for (cc in motif_cols) out2[[cc]] <- sum(out[[cc]], na.rm = TRUE)
  out2
}

.get_consensus_row <- function(gene_id,
                               condition,
                               final_consensus_path = "tables/final_consensus_table.csv") {
  tab <- .safe_read_csv(final_consensus_path)
  if (is.null(tab)) return(NULL)
  if (!all(c("AGI", "condition") %in% colnames(tab))) return(NULL)
  
  tab$AGI <- .clean_agi(tab$AGI)
  tab$condition <- as.character(tab$condition)
  
  gene_id <- .clean_agi(gene_id)
  out <- tab[tab$AGI == gene_id & tab$condition == condition, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)
  out[1, , drop = FALSE]
}

## "Defensible strength": min/median padj per method across SUPPORTING gears
## Supporting gears are taken from tables/support_by_gene_condition_gear.csv:
##  - robust: gears where support_votes_strict == TRUE
##  - candidate-only: gears where support_votes_sugg   == TRUE
.get_rhythmicity_strength <- function(gene_id,
                                      condition,
                                      final_consensus_path = "tables/final_consensus_table.csv",
                                      support_by_gear_path = "tables/support_by_gene_condition_gear.csv",
                                      method_results_path  = "tables/method_results_long.csv",
                                      vote_methods = c("RAIN", "JTK", "COSINOR")) {
  gene_id <- .clean_agi(gene_id)
  
  # Required tables
  cons <- .get_consensus_row(gene_id, condition, final_consensus_path = final_consensus_path)
  if (is.null(cons)) return(NULL)
  
  sup <- .safe_read_csv(support_by_gear_path)
  meth <- .safe_read_csv(method_results_path)
  
  if (is.null(sup) || is.null(meth)) {
    # If these are gitignored, don't crash; just return NULL
    return(NULL)
  }
  
  # Basic schema checks
  need_sup <- c("AGI","condition","gear_id","support_votes_strict","support_votes_sugg")
  need_meth <- c("AGI","condition","gear_id","method","padj","rhythmic_strict","rhythmic_suggestive")
  if (!all(need_sup %in% colnames(sup))) return(NULL)
  if (!all(need_meth %in% colnames(meth))) return(NULL)
  
  sup$AGI <- .clean_agi(sup$AGI)
  sup$condition <- as.character(sup$condition)
  
  meth$AGI <- .clean_agi(meth$AGI)
  meth$condition <- as.character(meth$condition)
  meth$method <- toupper(as.character(meth$method))
  
  # Determine classification for choosing supporting gears
  is_robust <- isTRUE(cons$consensus_rhythmic_phase) || isTRUE(cons$consensus_rhythmic)
  is_cand   <- isTRUE(cons$consensus_candidate)
  
  tier <- if (is_robust) "ROBUST" else if (is_cand) "CANDIDATE" else "NOT_OSCILLATOR"
  
  sup_g <- sup[sup$AGI == gene_id & sup$condition == condition, , drop = FALSE]
  if (nrow(sup_g) == 0) return(data.frame(
    AGI = gene_id, condition = condition, tier = tier,
    method = character(0), n_gears = integer(0),
    min_padj = numeric(0), median_padj = numeric(0),
    n_strict = integer(0), n_suggestive = integer(0),
    stringsAsFactors = FALSE
  ))
  
  # Choose supporting gears depending on tier
  if (tier == "ROBUST") {
    gears_keep <- unique(sup_g$gear_id[isTRUE(sup_g$support_votes_strict) | sup_g$support_votes_strict == TRUE])
  } else if (tier == "CANDIDATE") {
    gears_keep <- unique(sup_g$gear_id[isTRUE(sup_g$support_votes_sugg) | sup_g$support_votes_sugg == TRUE])
  } else {
    gears_keep <- integer(0)
  }
  
  # Pull method-level padj for those gears
  meth_g <- meth[
    meth$AGI == gene_id &
      meth$condition == condition &
      meth$gear_id %in% gears_keep &
      meth$method %in% vote_methods &
      is.finite(meth$padj),
    ,
    drop = FALSE
  ]
  
  # If no supporting gears (or no finite padj), return structured empty result
  if (length(gears_keep) == 0 || nrow(meth_g) == 0) {
    return(data.frame(
      AGI = gene_id, condition = condition, tier = tier,
      method = vote_methods,
      n_gears = 0L,
      min_padj = NA_real_,
      median_padj = NA_real_,
      n_strict = 0L,
      n_suggestive = 0L,
      stringsAsFactors = FALSE
    ))
  }
  
  # Summarise padj by method across supporting gears
  out <- lapply(vote_methods, function(mm) {
    x <- meth_g[meth_g$method == mm, , drop = FALSE]
    if (nrow(x) == 0) {
      return(data.frame(
        AGI = gene_id, condition = condition, tier = tier,
        method = mm,
        n_gears = 0L,
        min_padj = NA_real_,
        median_padj = NA_real_,
        n_strict = 0L,
        n_suggestive = 0L,
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      AGI = gene_id,
      condition = condition,
      tier = tier,
      method = mm,
      n_gears = length(unique(x$gear_id)),
      min_padj = min(x$padj, na.rm = TRUE),
      median_padj = median(x$padj, na.rm = TRUE),
      n_strict = sum(x$rhythmic_strict %in% TRUE, na.rm = TRUE),
      n_suggestive = sum(x$rhythmic_suggestive %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  
  # Optional "combined" summaries (not inferential, but transparent)
  out$neglog10_median <- suppressWarnings(-log10(out$median_padj))
  out
}

## ----------------------------
## 3) Raw plotting (your original behavior)
## ----------------------------

plot_gene_condition_raw <- function(gene_id,
                                    expr_mat,
                                    annot,
                                    condition_label = "condition") {
  .validate_inputs(gene_id, expr_mat, annot)
  
  gene_long    <- .build_gene_long(gene_id, expr_mat, annot)
  gene_summary <- .summarise_by_gear_time(gene_long)
  
  # Overlay: each GEAR separate color
  p_all <- ggplot(gene_summary,
                  aes(x = time, y = expression_mean,
                      color = factor(gear_id), group = gear_id)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = paste0(gene_id, " • ", condition_label, " across datasets"),
      x = "Time",
      y = "MBQN expression (mean per time)",
      color = "GEAR ID"
    ) +
    theme_bw(base_size = 14)
  
  # Facets: one panel per gear
  p_facets <- ggplot(gene_summary, aes(x = time, y = expression_mean, group = 1)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(size = 2, color = "black") +
    facet_wrap(~ gear_id, scales = "free_x") +
    labs(
      title = paste0(gene_id, " • ", condition_label, " per dataset"),
      x = "Time",
      y = "MBQN expression (mean per time)"
    ) +
    theme_bw(base_size = 14)
  
  list(
    data_raw_long    = gene_long,
    data_raw_summary = gene_summary,
    plot_all         = p_all,
    plot_panels      = p_facets
  )
}

## ----------------------------
## 4) Normalized 0..1 overlay (gray gears + black mean)
## ----------------------------

plot_gene_condition_norm01 <- function(gene_id,
                                       expr_mat,
                                       annot,
                                       condition_label = "condition",
                                       min_range = 0.1,
                                       require_shared_times = FALSE,
                                       show_points = TRUE) {
  .validate_inputs(gene_id, expr_mat, annot)
  
  gene_long <- .build_gene_long(gene_id, expr_mat, annot)
  gene_avg  <- gene_long %>%
    group_by(gear_id, time) %>%
    summarise(expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
    arrange(gear_id, time)
  
  # Min–max normalize within each gear (robust to all-NA gears)
  df_norm <- gene_avg %>%
    group_by(gear_id) %>%
    mutate(
      n_ok = sum(is.finite(expr)),
      min_expr = ifelse(n_ok > 0, min(expr[is.finite(expr)], na.rm = TRUE), NA_real_),
      max_expr = ifelse(n_ok > 0, max(expr[is.finite(expr)], na.rm = TRUE), NA_real_),
      range    = max_expr - min_expr,
      expr_norm01 = ifelse(is.finite(range) & range >= min_range,
                           (expr - min_expr) / range,
                           NA_real_)
    ) %>%
    ungroup() %>%
    filter(is.finite(expr_norm01))
  
  if (nrow(df_norm) == 0) {
    stop("No usable series after normalization. Try lowering min_range (currently ", min_range, ").")
  }
  
  # Optional: only keep timepoints shared across all gears kept
  if (require_shared_times) {
    gears_kept <- sort(unique(df_norm$gear_id))
    shared_times <- df_norm %>%
      count(gear_id, time) %>%
      count(time) %>%
      filter(n == length(gears_kept)) %>%
      pull(time)
    
    df_norm <- df_norm %>% filter(time %in% shared_times)
  }
  
  df_mean <- df_norm %>%
    group_by(time) %>%
    summarise(
      mean_norm01 = mean(expr_norm01, na.rm = TRUE),
      n_gears     = n_distinct(gear_id),
      .groups = "drop"
    ) %>%
    arrange(time)
  
  p <- ggplot(df_norm, aes(x = time, y = expr_norm01, group = factor(gear_id))) +
    geom_line(color = "grey70", linewidth = 0.9, alpha = 0.9) +
    { if (show_points) geom_point(color = "grey45", size = 1.8, alpha = 0.7) } +
    geom_line(
      data = df_mean,
      mapping = aes(x = time, y = mean_norm01),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1.3
    ) +
    geom_point(
      data = df_mean,
      mapping = aes(x = time, y = mean_norm01),
      inherit.aes = FALSE,
      color = "black",
      size = 2.2
    ) +
    labs(
      title = paste0(gene_id, " • ", condition_label, " (normalized 0–1)"),
      subtitle = if (require_shared_times)
        "Gray = per GEAR (min–max); Black = mean across gears (shared timepoints only)"
      else
        "Gray = per GEAR (min–max); Black = mean across gears (timepoints pooled)",
      x = "Time",
      y = "Normalized expression (0=min, 1=max)"
    ) +
    theme_bw(base_size = 14)
  
  list(
    data_norm_pergear = df_norm,
    data_norm_mean    = df_mean,
    plot_norm_overlay = p
  )
}

## ----------------------------
## 5) One wrapper: choose mode (+ optional attachments)
## ----------------------------

plot_gene_by_condition <- function(gene_id,
                                   condition = c("12L12D", "12L12D_LL", "16L8D", "8L16D"),
                                   mode = c("raw", "norm01"),
                                   min_range = 0.1,
                                   require_shared_times = FALSE,
                                   attach_motifs = TRUE,
                                   attach_strength = TRUE,
                                   print_summaries = TRUE) {
  cond <- match.arg(condition)
  mode <- match.arg(mode)
  
  gene_id_clean <- .clean_agi(gene_id)
  dat <- get_condition_data(cond)
  
  res <- if (mode == "raw") {
    plot_gene_condition_raw(
      gene_id = gene_id_clean,
      expr_mat = dat$expr,
      annot = dat$annot,
      condition_label = dat$label
    )
  } else {
    plot_gene_condition_norm01(
      gene_id = gene_id_clean,
      expr_mat = dat$expr,
      annot = dat$annot,
      condition_label = dat$label,
      min_range = min_range,
      require_shared_times = require_shared_times
    )
  }
  
  # Attach motif counts (this condition only)
  if (attach_motifs) {
    res$motif_summary <- .get_motif_summary(gene_id_clean, cond)
  }
  
  # Attach rhythmicity strength (padj summaries across supporting gears)
  if (attach_strength) {
    res$rhythmicity_strength <- .get_rhythmicity_strength(gene_id_clean, cond)
  }
  
  # Optional: print summaries
  if (print_summaries) {
    cat("\n============================\n")
    cat("Gene:", gene_id_clean, "\n")
    cat("Condition:", cond, "\n")
    cat("Mode:", mode, "\n")
    
    if (!is.null(res$rhythmicity_strength)) {
      cat("\n[Rhythmicity strength: padj across SUPPORTING gears]\n")
      print(res$rhythmicity_strength, row.names = FALSE)
    } else {
      cat("\n[Rhythmicity strength]\n")
      cat("  (Not available. Expected files:\n")
      cat("   - tables/final_consensus_table.csv\n")
      cat("   - tables/support_by_gene_condition_gear.csv\n")
      cat("   - tables/method_results_long.csv )\n")
    }
    
    if (!is.null(res$motif_summary)) {
      cat("\n[Motif counts in promoter (-1000..-1; both strands)]\n")
      print(res$motif_summary, row.names = FALSE)
    } else {
      cat("\n[Motif counts]\n")
      cat("  (Not available. Expected file:\n")
      cat("   - tables/motif_scan_per_photoperiod/motif_counts_", cond, ".csv )\n", sep = "")
    }
    cat("============================\n\n")
  }
  
  res
}

## ----------------------------
## 6) Example usage (interactive only)
## ----------------------------

if (interactive()) {
  gene_id   <- "AT3G57430"   # example gene
  condition <- "16L8D"
  
  # Raw plots (overlay + facets)
  res <- plot_gene_by_condition(gene_id, condition, mode = "raw")
  print(res$plot_all)
  print(res$plot_panels)
  
  # Normalized overlay (gray gears + black mean)
  resn <- plot_gene_by_condition(gene_id, condition, mode = "norm01",
                                 min_range = 0.1, require_shared_times = FALSE)
  print(resn$plot_norm_overlay)
}
