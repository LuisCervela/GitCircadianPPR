## 05_rain_rhythmicity_all.R
## Run RAIN on all conditions (12L12D, 12L12D→LL, 16L8D, 8L16D)
## and count rhythmic vs non-rhythmic PPR/clock genes.

library(rain)
library(dplyr)

## ---------- 1. Helper: run RAIN for a single GEAR ID ----------

run_rain_for_dataset <- function(gear,
                                 expr_mat,
                                 annot,
                                 period = 24,
                                 adj_method = "BH") {
  # Subset annotation for this dataset
  ann_g <- annot[annot$gear_id == gear, ]
  if (nrow(ann_g) == 0) {
    stop(paste("No samples found for GEAR", gear))
  }
  
  # Order by time
  ann_g <- ann_g[order(ann_g$time), ]
  times_full <- ann_g$time
  cols_full  <- ann_g$column_name
  
  # ---- 1. Average replicates at the same time ----
  times_unique <- sort(unique(times_full))
  
  X_avg_list <- lapply(times_unique, function(t) {
    cols_t <- cols_full[times_full == t]
    mat_t  <- expr_mat[, cols_t, drop = FALSE]   # genes x replicates (or 1)
    rowMeans(mat_t, na.rm = TRUE)               # works for 1 or >1 cols
  })
  
  X_avg <- do.call(cbind, X_avg_list)  # genes x timepoints
  rownames(X_avg) <- rownames(expr_mat)
  colnames(X_avg) <- paste0("t", times_unique)
  
  # ---- 2. Check time intervals using unique times ----
  dt <- diff(times_unique)
  dt_unique <- unique(dt[is.finite(dt)])
  
  if (length(times_unique) < 3) {
    message("GEAR ", gear, " has fewer than 3 unique timepoints. Skipping RAIN.")
    return(data.frame(
      AGI      = rownames(expr_mat),
      gear_id  = gear,
      p_value  = NA_real_,
      phase    = NA_real_,
      padj     = NA_real_,
      rhythmic = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  if (length(dt_unique) != 1) {
    message("GEAR ", gear,
            " has non-uniform time intervals (after averaging replicates): ",
            paste(dt_unique, collapse = ", "),
            ". Skipping RAIN for this dataset.")
    return(data.frame(
      AGI      = rownames(expr_mat),
      gear_id  = gear,
      p_value  = NA_real_,
      phase    = NA_real_,
      padj     = NA_real_,
      rhythmic = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  deltat <- dt_unique[1]
  
  # ---- 3. Prepare matrix for RAIN: rows = timepoints, cols = genes ----
  X_rain <- t(X_avg)   # timepoints x genes
  
  cat("Running RAIN for GEAR", gear,
      "| timepoints:", length(times_unique),
      "| deltat:", deltat, "\n")
  
  # ---- 4. Call rain() safely ----
  res <- tryCatch(
    rain(
      x      = X_rain,
      deltat = deltat,
      period = period,
      method = "independent",
      na.rm  = TRUE
    ),
    error = function(e) {
      message("RAIN failed for GEAR ", gear, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  
  if (is.null(res)) {
    return(data.frame(
      AGI      = rownames(expr_mat),
      gear_id  = gear,
      p_value  = NA_real_,
      phase    = NA_real_,
      padj     = NA_real_,
      rhythmic = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  pvals  <- as.numeric(res$pVal)
  phases <- as.numeric(res$phase)
  
  if (length(pvals) != nrow(expr_mat)) {
    warning("Length of pvals != number of genes for GEAR ", gear,
            " (", length(pvals), " vs ", nrow(expr_mat), ")")
  }
  
  df <- data.frame(
    AGI      = rownames(expr_mat),
    gear_id  = gear,
    p_value  = pvals,
    phase    = phases,
    stringsAsFactors = FALSE
  )
  
  df$padj     <- p.adjust(df$p_value, method = adj_method)
  df$rhythmic <- df$padj < 0.05
  
  df
}

## ---------- 2. Summarize rhythmicity per condition ----------

summarize_rain_condition <- function(rain_df, all_genes) {
  # remove rows without p_value
  rain_df <- rain_df %>%
    filter(!is.na(p_value))
  
  summary_df <- rain_df %>%
    group_by(AGI) %>%
    summarise(
      n_datasets_tested   = n(),
      n_rhythmic_datasets = sum(rhythmic, na.rm = TRUE),
      rhythmic_any        = n_rhythmic_datasets >= 1,
      rhythmic_in_2plus   = n_rhythmic_datasets >= 2,
      rhythmic_in_3plus   = n_rhythmic_datasets >= 3,
      rhythmic_in_all     = n_rhythmic_datasets == n_datasets_tested,
      .groups = "drop"
    )
  
  # ensure all genes present
  missing_genes <- setdiff(all_genes, summary_df$AGI)
  if (length(missing_genes) > 0) {
    add_df <- data.frame(
      AGI                  = missing_genes,
      n_datasets_tested    = 0L,
      n_rhythmic_datasets  = 0L,
      rhythmic_any         = FALSE,
      rhythmic_in_2plus    = FALSE,
      rhythmic_in_3plus    = FALSE,
      rhythmic_in_all      = FALSE,
      stringsAsFactors     = FALSE
    )
    summary_df <- bind_rows(summary_df, add_df)
  }
  
  summary_df %>% arrange(AGI)
}

## ---------- 3. Helper: count TRUE/FALSE for a logical flag ----------

count_rhythmic_flags <- function(summary_df, flag_col = "rhythmic_in_3") {
  stopifnot(flag_col %in% colnames(summary_df))
  
  x <- summary_df[[flag_col]]
  
  data.frame(
    flag      = flag_col,
    n_TRUE    = sum(x == TRUE,  na.rm = TRUE),
    n_FALSE   = sum(x == FALSE, na.rm = TRUE),
    n_NA      = sum(is.na(x)),
    total     = length(x),
    stringsAsFactors = FALSE
  )
}

## ---------- 4. Check that loaded objects exist ----------

if (!exists("expr_12L12D")    || !exists("annot_12L12D")   ||
    !exists("expr_12L12D_LL") || !exists("annot_12L12D_LL")||
    !exists("expr_16L8D")     || !exists("annot_16L8D")    ||
    !exists("expr_8L16D")     || !exists("annot_8L16D")) {
  stop("Expression/annotation objects not found. Run 03_load_PPRclock_all_conditions.R first.")
}

## ---------- 5. 12L12D (LD) ----------

ld_gears <- sort(unique(annot_12L12D$gear_id))
cat("12L12D GEAR IDs:", paste(ld_gears, collapse = ", "), "\n")

rain_ld_list <- lapply(ld_gears, function(g) {
  run_rain_for_dataset(
    gear     = g,
    expr_mat = expr_12L12D,
    annot    = annot_12L12D,
    period   = 24
  )
})
rain_ld_results <- bind_rows(rain_ld_list)

ld_genes <- rownames(expr_12L12D)
ld_rain_summary <- summarize_rain_condition(
  rain_df   = rain_ld_results,
  all_genes = ld_genes
)

ld_counts <- count_rhythmic_flags(ld_rain_summary, "rhythmic_in_3plus")
ld_counts$condition <- "12L12D"
print(ld_counts)

## ---------- 6. 12L12D → LL ----------

ll_gears <- sort(unique(annot_12L12D_LL$gear_id))
cat("12L12D→LL GEAR IDs:", paste(ll_gears, collapse = ", "), "\n")

rain_ll_list <- lapply(ll_gears, function(g) {
  run_rain_for_dataset(
    gear     = g,
    expr_mat = expr_12L12D_LL,
    annot    = annot_12L12D_LL,
    period   = 24
  )
})
rain_ll_results <- bind_rows(rain_ll_list)

ll_genes <- rownames(expr_12L12D_LL)
ll_rain_summary <- summarize_rain_condition(
  rain_df   = rain_ll_results,
  all_genes = ll_genes
)

ll_counts <- count_rhythmic_flags(ll_rain_summary, "rhythmic_in_3plus")
ll_counts$condition <- "12L12D→LL"
print(ll_counts)

## ---------- 7. 16L8D (Long Day) ----------

ldlong_gears <- sort(unique(annot_16L8D$gear_id))
cat("16L8D GEAR IDs:", paste(ldlong_gears, collapse = ", "), "\n")

rain_16L8D_list <- lapply(ldlong_gears, function(g) {
  run_rain_for_dataset(
    gear     = g,
    expr_mat = expr_16L8D,
    annot    = annot_16L8D,
    period   = 24
  )
})
rain_16L8D_results <- bind_rows(rain_16L8D_list)

genes_16L8D <- rownames(expr_16L8D)
rain_16L8D_summary <- summarize_rain_condition(
  rain_df   = rain_16L8D_results,
  all_genes = genes_16L8D
)

counts_16L8D <- count_rhythmic_flags(rain_16L8D_summary, "rhythmic_in_3plus")
counts_16L8D$condition <- "16L8D"
print(counts_16L8D)

## ---------- 8. 8L16D (Short Day) ----------

sd_gears <- sort(unique(annot_8L16D$gear_id))
cat("8L16D GEAR IDs:", paste(sd_gears, collapse = ", "), "\n")

rain_8L16D_list <- lapply(sd_gears, function(g) {
  run_rain_for_dataset(
    gear     = g,
    expr_mat = expr_8L16D,
    annot    = annot_8L16D,
    period   = 24
  )
})
rain_8L16D_results <- bind_rows(rain_8L16D_list)

genes_8L16D <- rownames(expr_8L16D)
rain_8L16D_summary <- summarize_rain_condition(
  rain_df   = rain_8L16D_results,
  all_genes = genes_8L16D
)

counts_8L16D <- count_rhythmic_flags(rain_8L16D_summary, "rhythmic_in_3plus")
counts_8L16D$condition <- "8L16D"
print(counts_8L16D)

## ---------- 9. Combined counts table (optional) ----------

rain_counts_all <- bind_rows(
  ld_counts,
  ll_counts,
  counts_16L8D,
  counts_8L16D
) %>%
  select(condition, flag, n_TRUE, n_FALSE, n_NA, total)

print(rain_counts_all)






