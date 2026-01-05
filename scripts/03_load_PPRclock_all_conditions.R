## ============================================================
## 03_load_PPRclock_all_conditions.R
##
## Purpose
##   Load the PPR+clock subset tables (Script 02 outputs) for all photoperiod
##   conditions and create, in the global environment:
##
##     annot_12L12D,    expr_12L12D
##     annot_12L12D_LL, expr_12L12D_LL
##     annot_16L8D,     expr_16L8D
##     annot_8L16D,     expr_8L16D
##
## Each loaded CSV is assumed to have:
##   - first 3 rows = metadata:
##       row 1: GEAR IDs (per sample column)
##       row 2: time (numeric, per sample column)
##       row 3: accession (string, per sample column)
##   - remaining rows = expression matrix (genes x samples)
##
## Inputs (from Script 02)
##   - data_clean/MBQN_12L12D_PPRclock.csv
##   - data_clean/MBQN_12L12D_to_LL_PPRclock.csv
##   - data_clean/MBQN_16L8D_PPRclock.csv
##   - data_clean/MBQN_8L16D_PPRclock.csv
##
## Notes
##   - check.names=FALSE is used to preserve exact column names.
##   - Expression values are coerced to numeric.
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

## ----------------------------
## Helper: load one condition
## ----------------------------
load_condition_PPRclock <- function(path_csv) {
  if (!file.exists(path_csv)) {
    stop("Missing input file: ", path_csv)
  }
  
  dat <- fread(path_csv, data.table = FALSE, check.names = FALSE)
  
  if (nrow(dat) < 4) {
    stop("File too small; expected >= 4 rows (3 metadata + expression): ", path_csv)
  }
  
  ## 1) Extract metadata rows
  gear_row <- as.character(dat[1, , drop = TRUE])
  time_row <- as.character(dat[2, , drop = TRUE])
  acc_row  <- as.character(dat[3, , drop = TRUE])
  
  cn <- colnames(dat)
  if (length(cn) < 2) stop("Expected at least 2 columns in: ", path_csv)
  
  sample_cols <- cn[-1]  # skip first column ("GEAR ID" / gene id column)
  
  annot <- data.frame(
    column_name = sample_cols,
    gear_id     = suppressWarnings(as.numeric(sub("\\..*$", "", gear_row[-1]))),
    time        = suppressWarnings(as.numeric(time_row[-1])),
    accession   = acc_row[-1],
    stringsAsFactors = FALSE
  )
  
  ## 2) Expression matrix
  dat_expr <- dat[-c(1:3), , drop = FALSE]
  colnames(dat_expr)[1] <- "AGI"
  
  agi_vec <- as.character(dat_expr$AGI)
  
  # Align columns explicitly to annotation column names
  expr_df <- dat_expr[, annot$column_name, drop = FALSE]
  
  # Coerce to numeric safely
  expr_df[] <- lapply(expr_df, function(x) suppressWarnings(as.numeric(as.character(x))))
  
  expr_mat <- as.matrix(expr_df)
  rownames(expr_mat) <- agi_vec
  
  list(annot = annot, expr = expr_mat)
}

## ----------------------------
## Load all conditions
## ----------------------------
paths <- list(
  `12L12D`    = "data_clean/MBQN_12L12D_PPRclock.csv",
  `12L12D_LL` = "data_clean/MBQN_12L12D_to_LL_PPRclock.csv",
  `16L8D`     = "data_clean/MBQN_16L8D_PPRclock.csv",
  `8L16D`     = "data_clean/MBQN_8L16D_PPRclock.csv"
)

## 12L12D (LD)
ld_res <- load_condition_PPRclock(paths[["12L12D"]])
annot_12L12D <- ld_res$annot
expr_12L12D  <- ld_res$expr
cat("[12L12D] expr_12L12D:", paste(dim(expr_12L12D), collapse = " x "), "\n")

## 12L12D → LL
ldll_res <- load_condition_PPRclock(paths[["12L12D_LL"]])
annot_12L12D_LL <- ldll_res$annot
expr_12L12D_LL  <- ldll_res$expr
cat("[12L12D_LL] expr_12L12D_LL:", paste(dim(expr_12L12D_LL), collapse = " x "), "\n")

## 16L8D (Long Day)
ldlong_res <- load_condition_PPRclock(paths[["16L8D"]])
annot_16L8D <- ldlong_res$annot
expr_16L8D  <- ldlong_res$expr
cat("[16L8D] expr_16L8D:", paste(dim(expr_16L8D), collapse = " x "), "\n")

## 8L16D (Short Day)
sd_res <- load_condition_PPRclock(paths[["8L16D"]])
annot_8L16D <- sd_res$annot
expr_8L16D  <- sd_res$expr
cat("[8L16D] expr_8L16D:", paste(dim(expr_8L16D), collapse = " x "), "\n")

## ----------------------------
## Optional sanity check
## ----------------------------
test_gene <- "AT5G61380"  # TOC1

cat("Is", test_gene, "in 12L12D?     ", test_gene %in% rownames(expr_12L12D), "\n")
cat("Is", test_gene, "in 12L12D_LL?  ", test_gene %in% rownames(expr_12L12D_LL), "\n")
cat("Is", test_gene, "in 16L8D?      ", test_gene %in% rownames(expr_16L8D), "\n")
cat("Is", test_gene, "in 8L16D?      ", test_gene %in% rownames(expr_8L16D), "\n")