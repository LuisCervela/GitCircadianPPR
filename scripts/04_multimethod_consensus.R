## ============================================================
## 04_multimethod_consensus.R
##
## Purpose
##   Multi-method rhythmicity calling across multiple independent datasets (GEARs),
##   within each photoperiod condition, followed by consensus calling across GEARs.
##
##   Methods (per GEAR dataset, after averaging replicates at identical timepoints):
##     - RAIN (rain package)
##     - COSINOR (linear regression with sin/cos terms; scan period grid + Bonferroni)
##     - JTK (MetaCycle::meta2d)
##     - ARS (MetaCycle::meta2d)  [auxiliary; not used in voting by default]
##
##   Consensus rule (default):
##     support_votes_strict = RAIN_strict AND (JTK_strict OR COSINOR_strict)
##     support_votes_sugg   = RAIN_sugg   AND (JTK_sugg   OR COSINOR_sugg)
##
##   Replication requirement:
##     - compute eligible GEAR datasets per condition using time coverage criteria
##     - require >= gears_required supporting GEARs (adaptive up to min_gears_cap)
##
## Inputs (required in memory; run Script 03 first)
##   expr_12L12D,    annot_12L12D
##   expr_12L12D_LL, annot_12L12D_LL
##   expr_16L8D,     annot_16L8D
##   expr_8L16D,     annot_8L16D
##
## Outputs (written to tables/)
##   - eligibility_by_condition_gear.csv
##   - method_results_long.csv
##   - support_by_gene_condition_gear.csv
##   - consensus_by_gene_condition.csv
##   - final_consensus_table.csv
##   - intermediate MetaCycle files in tables/metacycle_tmp/
##
## Notes / assumptions
##   - annot objects must contain: column_name, gear_id, time (numeric), accession
##   - expr matrices are genes x samples, with colnames matching annot$column_name
##   - timepoints may have replicates; these are averaged per timepoint per GEAR
##   - require_uniform_dt enforces a single sampling interval within each GEAR
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rain)
  library(MetaCycle)
  library(lomb) # kept for compatibility; not directly used here
  library(tidyr)
  library(openxlsx)
})

## ----------------------------
## 0) USER PARAMETERS
## ----------------------------
period <- 24

# Eligibility thresholds (applied to UNIQUE timepoints per GEAR after averaging replicates)
min_timepoints      <- 12
min_span_hours      <- 40
require_uniform_dt  <- TRUE

# Replication requirement cap (adaptive)
min_gears_cap <- 3

# Methods used for voting (ARS is auxiliary by default)
vote_methods <- c("RAIN", "COSINOR", "JTK")

# COSINOR scanning grid (keep narrow to avoid inflation)
cosinor_period_grid <- 22:26

# MetaCycle filters
metacycle_require_complete <- TRUE
metacycle_drop_constant    <- TRUE
metacycle_sd_eps           <- 1e-8

# Phase-consistency QC window (hours)
phase_window <- 4

# Output dirs
out_tables <- "tables"
out_tmp_mc <- file.path(out_tables, "metacycle_tmp")
dir.create(out_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(out_tmp_mc, showWarnings = FALSE, recursive = TRUE)

## ----------------------------
## 0.0) Safety: required objects exist (from Script 03)
## ----------------------------
stopifnot(exists("expr_12L12D"),    exists("annot_12L12D"))
stopifnot(exists("expr_12L12D_LL"), exists("annot_12L12D_LL"))
stopifnot(exists("expr_16L8D"),     exists("annot_16L8D"))
stopifnot(exists("expr_8L16D"),     exists("annot_8L16D"))

## ----------------------------
## 0.1) Helper functions
## ----------------------------

add_calls <- function(df, p_col = "p_value", adj_method = "BH") {
  df$padj <- NA_real_
  ok <- is.finite(df[[p_col]])
  df$padj[ok] <- p.adjust(df[[p_col]][ok], method = adj_method)
  df$rhythmic_strict     <- ifelse(is.na(df$padj), NA, df$padj < 0.05)
  df$rhythmic_suggestive <- ifelse(is.na(df$padj), NA, df$padj < 0.10)
  df
}

circ_mean_hours <- function(h, period = 24) {
  h <- h[is.finite(h)]
  if (length(h) == 0) return(NA_real_)
  ang <- 2*pi*h/period
  m <- atan2(mean(sin(ang)), mean(cos(ang)))
  if (!is.finite(m)) return(NA_real_)
  (m %% (2*pi)) * (period/(2*pi))
}

phase_dist_hours <- function(a, b, period = 24) {
  d <- (a - b + period/2) %% period - period/2
  abs(d)
}

detrend_linear <- function(time, y) {
  ok <- is.finite(time) & is.finite(y)
  if (sum(ok) < 4) return(y)
  fit <- stats::lm(y[ok] ~ time[ok])
  y2 <- y
  y2[ok] <- stats::residuals(fit) + mean(y[ok], na.rm = TRUE)
  y2
}

## ----------------------------
## 0.2) Condition object getter (requires Script 03 objects loaded)
## ----------------------------
get_condition_objects <- function(condition) {
  cond <- match.arg(condition, c("12L12D","12L12D_LL","16L8D","8L16D"))
  if (cond == "12L12D")    return(list(expr=expr_12L12D,    annot=annot_12L12D,    label="12L12D"))
  if (cond == "12L12D_LL") return(list(expr=expr_12L12D_LL, annot=annot_12L12D_LL, label="12L12D_LL"))
  if (cond == "16L8D")     return(list(expr=expr_16L8D,     annot=annot_16L8D,     label="16L8D"))
  if (cond == "8L16D")     return(list(expr=expr_8L16D,     annot=annot_8L16D,     label="8L16D"))
}

## ----------------------------
## 0.3) Average replicates at same time (per GEAR)
## ----------------------------
average_reps_by_time <- function(expr_mat, annot_g) {
  annot_g <- annot_g[order(annot_g$time), , drop = FALSE]
  times_full <- annot_g$time
  cols_full  <- annot_g$column_name
  times_unique <- sort(unique(times_full))
  
  X_avg_list <- lapply(times_unique, function(t) {
    cols_t <- cols_full[times_full == t]
    mat_t  <- expr_mat[, cols_t, drop = FALSE]
    rowMeans(mat_t, na.rm = TRUE)
  })
  
  X_avg <- do.call(cbind, X_avg_list)
  rownames(X_avg) <- rownames(expr_mat)
  colnames(X_avg) <- paste0("t", times_unique)
  
  dt <- if (length(times_unique) >= 2) diff(times_unique) else numeric(0)
  dt_unique <- unique(dt[is.finite(dt)])
  dt_val <- if (length(dt_unique) == 1) dt_unique[1] else NA_real_
  
  list(times_unique = times_unique, X_avg = X_avg, dt_unique = dt_unique, dt = dt_val)
}

## ============================================================
## 1) Eligibility table (per condition × GEAR)
## ============================================================
build_elig_table <- function(condition) {
  obj <- get_condition_objects(condition)
  expr <- obj$expr
  annot <- obj$annot
  gears <- sort(unique(annot$gear_id))
  
  out <- lapply(gears, function(g) {
    ann_g <- annot[annot$gear_id == g, , drop = FALSE]
    avg <- average_reps_by_time(expr, ann_g)
    times <- avg$times_unique
    span <- if (length(times) >= 2) (max(times) - min(times)) else 0
    dt_unique <- avg$dt_unique
    dt_ok <- if (!require_uniform_dt) TRUE else (length(dt_unique) == 1)
    
    eligible <- (length(times) >= min_timepoints) && (span >= min_span_hours) && dt_ok
    
    data.frame(
      condition = obj$label,
      gear_id = g,
      eligible = eligible,
      n_unique_time = length(times),
      t_min = if (length(times) > 0) min(times) else NA_real_,
      t_max = if (length(times) > 0) max(times) else NA_real_,
      span = span,
      dt = ifelse(length(dt_unique) == 1, dt_unique[1], NA_real_),
      effective_span = span + ifelse(length(dt_unique) == 1, dt_unique[1], NA_real_),
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(out)
}

conditions <- c("12L12D","12L12D_LL","16L8D","8L16D")

elig_table <- dplyr::bind_rows(lapply(conditions, build_elig_table))
write.csv(elig_table, file.path(out_tables, "eligibility_by_condition_gear.csv"), row.names = FALSE)

gears_required_by_condition <- elig_table %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    n_eligible_gears_condition = sum(eligible == TRUE, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    gears_required = dplyr::case_when(
      n_eligible_gears_condition <= 1 ~ 1L,
      n_eligible_gears_condition <= 3 ~ 2L,
      TRUE ~ as.integer(min_gears_cap)
    ),
    gears_required = pmin(gears_required, n_eligible_gears_condition)
  )

## ============================================================
## 2) RAIN per dataset
## ============================================================
run_rain_for_dataset <- function(gear, expr_mat, annot, period = 24, adj_method = "BH") {
  ann_g <- annot[annot$gear_id == gear, , drop = FALSE]
  if (nrow(ann_g) == 0) stop("No samples for GEAR ", gear)
  
  avg <- average_reps_by_time(expr_mat, ann_g)
  times_unique <- avg$times_unique
  X_avg <- avg$X_avg
  dt_unique <- avg$dt_unique
  
  if (length(times_unique) < 3 || (require_uniform_dt && length(dt_unique) != 1)) {
    df <- data.frame(
      AGI = rownames(expr_mat), gear_id = gear,
      p_value = NA_real_, phase = NA_real_,
      n_timepoints = length(times_unique), deltat = NA_real_,
      stringsAsFactors = FALSE
    )
    return(add_calls(df, adj_method = adj_method))
  }
  
  deltat <- dt_unique[1]
  X_rain <- t(X_avg)
  
  res <- tryCatch(
    rain::rain(x = X_rain, deltat = deltat, period = period, method = "independent", na.rm = TRUE),
    error = function(e) NULL
  )
  
  if (is.null(res)) {
    df <- data.frame(
      AGI = rownames(expr_mat), gear_id = gear,
      p_value = NA_real_, phase = NA_real_,
      n_timepoints = length(times_unique), deltat = deltat,
      stringsAsFactors = FALSE
    )
    return(add_calls(df, adj_method = adj_method))
  }
  
  df <- data.frame(
    AGI = rownames(expr_mat),
    gear_id = gear,
    p_value = as.numeric(res$pVal),
    phase = as.numeric(res$phase),
    n_timepoints = length(times_unique),
    deltat = deltat,
    stringsAsFactors = FALSE
  )
  add_calls(df, adj_method = adj_method)
}

run_rain_condition <- function(condition, period = 24) {
  obj <- get_condition_objects(condition)
  expr <- obj$expr
  annot <- obj$annot
  
  gears_ok <- elig_table %>%
    dplyr::filter(condition == obj$label, eligible == TRUE) %>%
    dplyr::pull(gear_id) %>% unique() %>% sort()
  
  out <- lapply(gears_ok, function(g) {
    run_rain_for_dataset(g, expr, annot, period = period) %>%
      dplyr::mutate(condition = obj$label, method = "RAIN") %>%
      dplyr::select(AGI, condition, gear_id, method,
                    p_value, padj, rhythmic_strict, rhythmic_suggestive,
                    phase, n_timepoints, deltat)
  })
  
  dplyr::bind_rows(out)
}

## ============================================================
## 3) COSINOR per dataset (detrended; scan period grid)
## ============================================================
cosinor_fit_p <- function(time, y, period = 24) {
  ok <- is.finite(time) & is.finite(y)
  time <- time[ok]; y <- y[ok]
  if (length(y) < 4) return(list(p = NA_real_, phase = NA_real_))
  
  omega <- 2*pi/period
  cos_t <- cos(omega*time)
  sin_t <- sin(omega*time)
  
  m0 <- stats::lm(y ~ 1)
  m1 <- stats::lm(y ~ cos_t + sin_t)
  a  <- stats::anova(m0, m1)
  p  <- suppressWarnings(as.numeric(a$`Pr(>F)`[2]))
  
  b <- stats::coef(m1)
  bc <- unname(b["cos_t"])
  bs <- unname(b["sin_t"])
  if (!is.finite(bc) || !is.finite(bs)) return(list(p = p, phase = NA_real_))
  
  phi <- atan2(bs, bc)
  phase <- ((2*pi - phi) %% (2*pi)) * (period/(2*pi))
  list(p = p, phase = phase)
}

run_cosinor_for_dataset <- function(gear, expr_mat, annot,
                                    period_grid = 22:26,
                                    adj_method = "BH") {
  ann_g <- annot[annot$gear_id == gear, , drop = FALSE]
  if (nrow(ann_g) == 0) stop("No samples for GEAR ", gear)
  
  avg <- average_reps_by_time(expr_mat, ann_g)
  times_unique <- avg$times_unique
  X_avg <- avg$X_avg
  
  genes <- rownames(expr_mat)
  p_raw <- rep(NA_real_, length(genes))
  phase_best <- rep(NA_real_, length(genes))
  
  m <- length(period_grid)
  
  for (i in seq_along(genes)) {
    y <- detrend_linear(times_unique, as.numeric(X_avg[i, ]))
    
    pv <- rep(NA_real_, m)
    ph <- rep(NA_real_, m)
    for (k in seq_along(period_grid)) {
      fit <- cosinor_fit_p(times_unique, y, period = period_grid[k])
      pv[k] <- fit$p
      ph[k] <- fit$phase
    }
    
    if (all(!is.finite(pv))) {
      p_raw[i] <- NA_real_
      phase_best[i] <- NA_real_
    } else {
      kbest <- which.min(pv)
      pmin <- pv[kbest]
      p_raw[i] <- min(1, pmin * m)  # Bonferroni over scanned periods
      phase_best[i] <- ph[kbest]
    }
  }
  
  df <- data.frame(
    AGI = genes,
    gear_id = gear,
    p_value = p_raw,
    phase = phase_best,
    n_timepoints = length(times_unique),
    deltat = if (length(times_unique) >= 2) median(diff(times_unique)) else NA_real_,
    stringsAsFactors = FALSE
  )
  
  add_calls(df, p_col = "p_value", adj_method = adj_method)
}

run_cosinor_condition <- function(condition, period = 24) {
  obj <- get_condition_objects(condition)
  expr <- obj$expr
  annot <- obj$annot
  
  gears_ok <- elig_table %>%
    dplyr::filter(condition == obj$label, eligible == TRUE) %>%
    dplyr::pull(gear_id) %>% unique() %>% sort()
  
  out <- lapply(gears_ok, function(g) {
    run_cosinor_for_dataset(g, expr, annot, period_grid = cosinor_period_grid) %>%
      dplyr::mutate(condition = obj$label, method = "COSINOR") %>%
      dplyr::select(AGI, condition, gear_id, method,
                    p_value, padj, rhythmic_strict, rhythmic_suggestive,
                    phase, n_timepoints, deltat)
  })
  
  dplyr::bind_rows(out)
}

## ============================================================
## 4) MetaCycle runner (JTK, ARS) with detrending
## ============================================================
filter_for_metacycle <- function(X,
                                 require_complete = TRUE,
                                 drop_constant = TRUE,
                                 sd_eps = 1e-8) {
  Xf <- X
  
  if (require_complete) {
    keep <- rowSums(!is.finite(Xf)) == 0
    Xf <- Xf[keep, , drop = FALSE]
  }
  
  if (drop_constant && nrow(Xf) > 0) {
    sds <- apply(Xf, 1, sd, na.rm = TRUE)
    keep <- is.finite(sds) & (sds > sd_eps)
    Xf <- Xf[keep, , drop = FALSE]
  }
  
  Xf
}

write_metacycle_infile <- function(X, times, path) {
  dat <- data.frame(CycID = rownames(X), X, check.names = FALSE)
  colnames(dat)[-1] <- as.character(times)
  write.table(dat, file = path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  invisible(path)
}

run_metacycle_onegear <- function(X, times, out_dir,
                                  method = c("JTK","ARS"),
                                  period = 24,
                                  quiet = TRUE) {
  method <- match.arg(method)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  infile <- file.path(out_dir, "metacycle_input.txt")
  write_metacycle_infile(X, times, infile)
  
  if (quiet) {
    invisible(capture.output({
      invisible(capture.output({
        MetaCycle::meta2d(
          infile = infile,
          filestyle = "txt",
          timepoints = times,
          cycMethod = method,
          minper = period,
          maxper = period,
          outdir = out_dir,
          outputFile = TRUE,
          outIntegration = "noIntegration"
        )
      }, type = "message"))
    }))
  } else {
    MetaCycle::meta2d(
      infile = infile,
      filestyle = "txt",
      timepoints = times,
      cycMethod = method,
      minper = period,
      maxper = period,
      outdir = out_dir,
      outputFile = TRUE,
      outIntegration = "noIntegration"
    )
  }
  
  list.files(out_dir, pattern = "\\.txt$", full.names = TRUE)
}

read_metacycle_result <- function(files, method = c("JTK","ARS")) {
  method <- match.arg(method)
  if (length(files) == 0) return(NULL)
  
  f_method <- files[grepl(method, basename(files), ignore.case = TRUE)]
  if (length(f_method) == 0) f_method <- files
  
  info <- file.info(f_method)
  f_main <- rownames(info)[which.max(info$size)]
  
  tab <- tryCatch(
    read.table(f_main, header = TRUE, sep = "\t",
               stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(tab) || nrow(tab) == 0) return(NULL)
  
  idcol <- if ("CycID" %in% colnames(tab)) "CycID" else colnames(tab)[1]
  tab$CycID <- tab[[idcol]]
  
  pick_first <- function(cands) {
    cands <- cands[cands %in% colnames(tab)]
    if (length(cands) == 0) return(NA_character_)
    cands[1]
  }
  
  if (method == "ARS") {
    pcol    <- pick_first(c("pvalue","P","p","Pvalue","P.Value"))
    padjcol <- pick_first(c("fdr_BH","BH.Q","ADJ.P","adj.p","padj"))
    phcol   <- pick_first(c("phase","Phase","PHASE"))
  } else {
    pcol    <- pick_first(c("pvalue","P","p","Pvalue","P.Value"))
    padjcol <- pick_first(c("ADJ.P","BH.Q","adj.p","padj"))
    phcol   <- pick_first(c("LAG","lag","Phase","phase"))
  }
  
  data.frame(
    AGI     = tab$CycID,
    p_value = if (!is.na(pcol))    suppressWarnings(as.numeric(tab[[pcol]]))     else NA_real_,
    padj    = if (!is.na(padjcol)) suppressWarnings(as.numeric(tab[[padjcol]])) else NA_real_,
    phase   = if (!is.na(phcol))   suppressWarnings(as.numeric(tab[[phcol]]))   else NA_real_,
    stringsAsFactors = FALSE
  )
}

run_metacycle_condition <- function(condition, method = c("JTK","ARS"),
                                    period = 24, adj_method = "BH") {
  method <- match.arg(method)
  
  obj <- get_condition_objects(condition)
  expr <- obj$expr
  annot <- obj$annot
  
  gears_ok <- elig_table %>%
    dplyr::filter(condition == obj$label, eligible == TRUE) %>%
    dplyr::pull(gear_id) %>% unique() %>% sort()
  
  out_all <- vector("list", length(gears_ok))
  
  for (i in seq_along(gears_ok)) {
    g <- gears_ok[i]
    
    ann_g <- annot[annot$gear_id == g, , drop = FALSE]
    avg <- average_reps_by_time(expr, ann_g)
    times <- avg$times_unique
    X <- avg$X_avg
    
    Xf <- filter_for_metacycle(
      X,
      require_complete = metacycle_require_complete,
      drop_constant = metacycle_drop_constant,
      sd_eps = metacycle_sd_eps
    )
    
    if (nrow(Xf) == 0) {
      df <- data.frame(
        AGI = rownames(expr), gear_id = g,
        p_value = NA_real_, phase = NA_real_,
        n_timepoints = length(times), deltat = avg$dt,
        stringsAsFactors = FALSE
      )
      df <- add_calls(df, adj_method = adj_method)
      
      out_all[[i]] <- df %>%
        dplyr::mutate(condition = obj$label, method = method) %>%
        dplyr::select(AGI, condition, gear_id, method,
                      p_value, padj, rhythmic_strict, rhythmic_suggestive,
                      phase, n_timepoints, deltat)
      next
    }
    
    # detrend each gene before MetaCycle
    Xf_det <- Xf
    for (ii in seq_len(nrow(Xf_det))) {
      Xf_det[ii, ] <- detrend_linear(times, as.numeric(Xf_det[ii, ]))
    }
    Xf <- Xf_det
    
    out_dir <- file.path(out_tmp_mc, obj$label, paste0("GEAR_", g), method)
    
    files <- tryCatch(
      run_metacycle_onegear(Xf, times, out_dir = out_dir, method = method, period = period, quiet = TRUE),
      error = function(e) character(0)
    )
    
    mc <- read_metacycle_result(files, method = method)
    
    if (is.null(mc)) {
      df <- data.frame(
        AGI = rownames(expr), gear_id = g,
        p_value = NA_real_, phase = NA_real_,
        n_timepoints = length(times), deltat = avg$dt,
        stringsAsFactors = FALSE
      )
      df <- add_calls(df, adj_method = adj_method)
      
      out_all[[i]] <- df %>%
        dplyr::mutate(condition = obj$label, method = method) %>%
        dplyr::select(AGI, condition, gear_id, method,
                      p_value, padj, rhythmic_strict, rhythmic_suggestive,
                      phase, n_timepoints, deltat)
      next
    }
    
    df <- data.frame(
      AGI = mc$AGI,
      gear_id = g,
      p_value = mc$p_value,
      padj    = mc$padj,  # MetaCycle adjusted (preferred if present)
      phase   = mc$phase,
      n_timepoints = length(times),
      deltat = avg$dt,
      stringsAsFactors = FALSE
    )
    
    # If MetaCycle didn't provide padj, fallback BH across genes once
    if (all(!is.finite(df$padj))) {
      ok <- is.finite(df$p_value)
      df$padj <- NA_real_
      df$padj[ok] <- p.adjust(df$p_value[ok], method = adj_method)
    }
    
    df$rhythmic_strict     <- ifelse(is.na(df$padj), NA, df$padj < 0.05)
    df$rhythmic_suggestive <- ifelse(is.na(df$padj), NA, df$padj < 0.10)
    
    out_all[[i]] <- df %>%
      dplyr::mutate(condition = obj$label, method = method) %>%
      dplyr::select(AGI, condition, gear_id, method,
                    p_value, padj, rhythmic_strict, rhythmic_suggestive,
                    phase, n_timepoints, deltat)
  }
  
  dplyr::bind_rows(out_all)
}

run_jtk_condition <- function(condition, period = 24) run_metacycle_condition(condition, method = "JTK", period = period)
run_ars_condition <- function(condition, period = 24) run_metacycle_condition(condition, method = "ARS", period = period)

## ============================================================
## 5) Run all conditions (progress counter) -> method_long2
## ============================================================
total_steps <- length(conditions) * 4
step <- 0

tick <- function(label, condition) {
  step <<- step + 1
  pct <- 100 * step / total_steps
  cat(sprintf("Progress: %3d/%3d (%.1f%%) | %s | %s\n",
              step, total_steps, pct, label, condition))
}

run_block <- function(method_label, fun) {
  res_list <- vector("list", length(conditions))
  for (i in seq_along(conditions)) {
    tick(method_label, conditions[i])
    res_list[[i]] <- fun(conditions[i], period = period)
  }
  dplyr::bind_rows(res_list)
}

method_long2 <- dplyr::bind_rows(
  run_block("RAIN",    run_rain_condition),
  run_block("COSINOR", run_cosinor_condition),
  run_block("JTK",     run_jtk_condition),
  run_block("ARS",     run_ars_condition)
)

write.csv(method_long2, file.path(out_tables, "method_results_long.csv"), row.names = FALSE)

## ============================================================
## 6) Support by gear (gene × condition × GEAR)
## ============================================================
method_long2 <- method_long2 %>%
  dplyr::left_join(
    elig_table %>% dplyr::select(condition, gear_id, eligible),
    by = c("condition","gear_id")
  )

support_by_gear <- method_long2 %>%
  dplyr::filter(eligible == TRUE, method %in% vote_methods, is.finite(padj)) %>%
  dplyr::group_by(AGI, condition, gear_id) %>%
  dplyr::summarise(
    # strict booleans
    rain_strict = any(method == "RAIN"    & rhythmic_strict == TRUE),
    jtk_strict  = any(method == "JTK"     & rhythmic_strict == TRUE),
    cos_strict  = any(method == "COSINOR" & rhythmic_strict == TRUE),
    
    # suggestive booleans
    rain_sugg = any(method == "RAIN"    & rhythmic_suggestive == TRUE),
    jtk_sugg  = any(method == "JTK"     & rhythmic_suggestive == TRUE),
    cos_sugg  = any(method == "COSINOR" & rhythmic_suggestive == TRUE),
    
    # phase vote = circular mean across rhythmic methods
    phase_vote_strict = circ_mean_hours(phase[rhythmic_strict == TRUE], period = period),
    phase_vote_sugg   = circ_mean_hours(phase[rhythmic_suggestive == TRUE], period = period),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    # voting rule: RAIN + one additional method (JTK or COSINOR)
    support_votes_strict = rain_strict & (jtk_strict | cos_strict),
    support_votes_sugg   = rain_sugg   & (jtk_sugg   | cos_sugg),
    
    # method counts per gear (useful for ranking downstream)
    n_methods_strict = as.integer(rain_strict) + as.integer(jtk_strict) + as.integer(cos_strict),
    n_methods_sugg   = as.integer(rain_sugg)   + as.integer(jtk_sugg)   + as.integer(cos_sugg)
  )

write.csv(support_by_gear, file.path(out_tables, "support_by_gene_condition_gear.csv"), row.names = FALSE)

## ============================================================
## 7) Consensus per gene × condition (replication across eligible GEARs)
## ============================================================
consensus_gene_condition <- support_by_gear %>%
  dplyr::group_by(AGI, condition) %>%
  dplyr::summarise(
    n_eligible_gears_tested      = dplyr::n_distinct(gear_id),
    n_gears_supported_strict     = sum(support_votes_strict, na.rm = TRUE),
    n_gears_supported_sugg       = sum(support_votes_sugg,   na.rm = TRUE),
    
    # method coverage (within eligible gears)
    n_gears_rain_strict = sum(rain_strict, na.rm = TRUE),
    n_gears_jtk_strict  = sum(jtk_strict,  na.rm = TRUE),
    n_gears_cos_strict  = sum(cos_strict,  na.rm = TRUE),
    
    n_gears_rain_sugg = sum(rain_sugg, na.rm = TRUE),
    n_gears_jtk_sugg  = sum(jtk_sugg,  na.rm = TRUE),
    n_gears_cos_sugg  = sum(cos_sugg,  na.rm = TRUE),
    
    # phase estimate across supporting gears (prefer strict, fallback suggestive)
    phase_est = {
      v1 <- phase_vote_strict[support_votes_strict & is.finite(phase_vote_strict)]
      if (length(v1) > 0) circ_mean_hours(v1, period) else {
        v2 <- phase_vote_sugg[support_votes_sugg & is.finite(phase_vote_sugg)]
        if (length(v2) > 0) circ_mean_hours(v2, period) else NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  dplyr::left_join(gears_required_by_condition, by = "condition") %>%
  dplyr::mutate(
    consensus_rhythmic  = n_gears_supported_strict >= gears_required,
    consensus_candidate = n_gears_supported_sugg   >= gears_required
  ) %>%
  dplyr::arrange(condition, AGI)

# Phase-consistency (strict-supported gears) within phase_window hours
phase_consistent_strict <- support_by_gear %>%
  dplyr::filter(support_votes_strict, is.finite(phase_vote_strict)) %>%
  dplyr::group_by(AGI, condition) %>%
  dplyr::summarise(
    phase_ref_strict = circ_mean_hours(phase_vote_strict, period = period),
    n_gears_supported_strict_phase4h =
      sum(phase_dist_hours(phase_vote_strict, phase_ref_strict, period = period) <= phase_window),
    .groups = "drop"
  )

consensus_gene_condition <- consensus_gene_condition %>%
  dplyr::left_join(phase_consistent_strict, by = c("AGI","condition")) %>%
  dplyr::mutate(
    n_gears_supported_strict_phase4h =
      ifelse(is.na(n_gears_supported_strict_phase4h), 0L, as.integer(n_gears_supported_strict_phase4h)),
    consensus_rhythmic_phase = n_gears_supported_strict_phase4h >= gears_required
  )

write.csv(consensus_gene_condition, file.path(out_tables, "consensus_by_gene_condition.csv"), row.names = FALSE)

## ============================================================
## 8) Final deliverable table
## ============================================================
final_consensus_table <- consensus_gene_condition %>%
  dplyr::select(
    condition, AGI,
    gears_required, n_eligible_gears_tested,
    n_gears_supported_strict, n_gears_supported_sugg,
    n_gears_supported_strict_phase4h,
    n_gears_rain_strict, n_gears_jtk_strict, n_gears_cos_strict,
    n_gears_rain_sugg,   n_gears_jtk_sugg,   n_gears_cos_sugg,
    phase_est, phase_ref_strict,
    consensus_rhythmic, consensus_rhythmic_phase,
    consensus_candidate
  ) %>%
  dplyr::arrange(condition, dplyr::desc(consensus_rhythmic), dplyr::desc(consensus_candidate), AGI)

write.csv(final_consensus_table, file.path(out_tables, "final_consensus_table.csv"), row.names = FALSE)

cat("Done. Final table written to:", file.path(out_tables, "final_consensus_table.csv"), "\n")

## ============================================================
## 9) Clock gene strength tables as a positive control
## ============================================================

clean_agi <- function(x) toupper(trimws(sub("\\..*$", "", x)))

# ---- inputs ----
clock <- readLines("data_clean/clock_only_genelist.txt") |> (\(x) x[nzchar(x)])()
clock <- clean_agi(clock)

final <- read.csv("tables/final_consensus_table.csv", stringsAsFactors = FALSE)
sup   <- read.csv("tables/support_by_gene_condition_gear.csv", stringsAsFactors = FALSE)
meth  <- read.csv("tables/method_results_long.csv", stringsAsFactors = FALSE)

final$AGI <- clean_agi(final$AGI)
sup$AGI   <- clean_agi(sup$AGI)
meth$AGI  <- clean_agi(meth$AGI)
meth$method <- toupper(as.character(meth$method))

vote_methods <- c("RAIN", "JTK", "COSINOR")
conditions <- c("12L12D", "12L12D_LL", "16L8D", "8L16D")

# ---- alias map (from your table) ----
alias_map <- tibble::tribble(
  ~AGI,        ~Alias,
  "AT5G64170", "LNK1",
  "AT5G61380", "TOC1",
  "AT3G46640", "LUX",
  "AT2G46830", "CCA1",
  "AT5G02810", "PRR7",
  "AT5G24470", "PRR5",
  "AT3G54500", "LNK2",
  "AT2G46790", "PRR9",
  "AT3G09600", "RVE8",
  "AT2G25930", "ELF3",
  "AT2G40080", "ELF4",
  "AT1G22770", "GI",
  "AT1G01060", "LHY"
) %>% mutate(AGI = clean_agi(AGI))

# ---- summariser (same semantics as Script 08 strength block) ----
strength_one <- function(gene_id, condition) {
  cons <- final %>% filter(AGI == gene_id, condition == condition)
  if (nrow(cons) == 0) return(NULL)
  cons <- cons[1, ]
  
  is_robust <- isTRUE(cons$consensus_rhythmic_phase) || isTRUE(cons$consensus_rhythmic)
  is_cand   <- isTRUE(cons$consensus_candidate)
  tier <- if (is_robust) "ROBUST" else if (is_cand) "CANDIDATE" else "NOT_OSCILLATOR"
  
  sup_g <- sup %>% filter(AGI == gene_id, condition == condition)
  
  gears_keep <- integer(0)
  if (nrow(sup_g) > 0) {
    if (tier == "ROBUST") {
      gears_keep <- unique(sup_g$gear_id[sup_g$support_votes_strict %in% TRUE])
    } else if (tier == "CANDIDATE") {
      gears_keep <- unique(sup_g$gear_id[sup_g$support_votes_sugg %in% TRUE])
    }
  }
  
  meth_g <- meth %>%
    filter(AGI == gene_id,
           condition == condition,
           gear_id %in% gears_keep,
           method %in% vote_methods,
           is.finite(padj))
  
  out <- lapply(vote_methods, function(mm) {
    x <- meth_g %>% filter(method == mm)
    if (nrow(x) == 0) {
      data.frame(
        AGI = gene_id, condition = condition, tier = tier,
        method = mm, n_gears = 0L,
        min_padj = NA_real_, median_padj = NA_real_,
        n_strict = 0L, n_suggestive = 0L
      )
    } else {
      data.frame(
        AGI = gene_id, condition = condition, tier = tier,
        method = mm,
        n_gears = n_distinct(x$gear_id),
        min_padj = min(x$padj, na.rm = TRUE),
        median_padj = median(x$padj, na.rm = TRUE),
        n_strict = sum(x$rhythmic_strict %in% TRUE, na.rm = TRUE),
        n_suggestive = sum(x$rhythmic_suggestive %in% TRUE, na.rm = TRUE)
      )
    }
  }) %>% bind_rows()
  
  out$neglog10_median <- suppressWarnings(-log10(out$median_padj))
  out
}

# ---- build tables ----
Clock_strength_long <- bind_rows(lapply(clock, function(g) {
  bind_rows(lapply(conditions, function(cc) strength_one(g, cc)))
})) %>%
  left_join(alias_map, by = "AGI") %>%
  relocate(Alias, .after = AGI) %>%
  arrange(condition, Alias, method)

Clock_strength_wide <- Clock_strength_long %>%
  select(AGI, Alias, condition, tier, method, median_padj) %>%
  pivot_wider(names_from = method, values_from = median_padj) %>%
  arrange(condition, Alias)

# ---- add phase estimates ----
phase_tab <- final %>%
  select(AGI, condition, phase_est) %>%
  mutate(phase_est = round(phase_est, 1))

Clock_strength_wide <- Clock_strength_wide %>%
  left_join(phase_tab, by = c("AGI","condition")) %>%
  relocate(phase_est, .after = tier)

# ---- export to Excel ----
wb <- createWorkbook()
addWorksheet(wb, "Clock13_strength_long")
addWorksheet(wb, "Clock13_strength_wide")

writeData(wb, "Clock13_strength_long", Clock_strength_long)
writeData(wb, "Clock13_strength_wide", Clock_strength_wide)

saveWorkbook(wb, "tables/Clock13_rhythmicity_strength.xlsx", overwrite = TRUE)
