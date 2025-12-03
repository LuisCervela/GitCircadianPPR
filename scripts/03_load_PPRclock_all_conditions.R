## 03_load_PPRclock_all_conditions.R
## Load PPR+clock tables for all conditions and build:
##  - annot_12L12D,    expr_12L12D
##  - annot_12L12D_LL, expr_12L12D_LL
##  - annot_16L8D,     expr_16L8D
##  - annot_8L16D,     expr_8L16D

library(data.table)
library(dplyr)

## ---------- Helper: load one condition ----------

load_condition_PPRclock <- function(path_csv) {
  dat <- fread(
    path_csv,
    data.table  = FALSE,
    check.names = FALSE
  )
  
  # 1) Extract metadata rows
  gear_row <- as.character(dat[1, ])
  time_row <- as.character(dat[2, ])
  acc_row  <- as.character(dat[3, ])
  
  cn <- colnames(dat)
  sample_cols <- cn[-1]  # skip first column ("GEAR ID" / "AGI")
  
  annot <- data.frame(
    column_name = sample_cols,
    gear_id     = suppressWarnings(as.numeric(sub("\\..*$", "", gear_row[-1]))),
    time        = suppressWarnings(as.numeric(time_row[-1])),
    accession   = acc_row[-1],
    stringsAsFactors = FALSE
  )
  
  # 2) Expression matrix
  dat_expr <- dat[-c(1:3), ]
  colnames(dat_expr)[1] <- "AGI"
  agi_vec <- dat_expr$AGI
  
  expr_df <- dat_expr[, annot$column_name, drop = FALSE]
  expr_df[] <- lapply(expr_df, function(x) as.numeric(as.character(x)))
  
  expr_mat <- as.matrix(expr_df)
  rownames(expr_mat) <- agi_vec
  
  list(
    annot = annot,
    expr  = expr_mat
  )
}

## ---------- 12L12D (LD) ----------

ld_res <- load_condition_PPRclock("data_clean/MBQN_12L12D_PPRclock.csv")
annot_12L12D <- ld_res$annot
expr_12L12D  <- ld_res$expr

cat("[12L12D] expr_12L12D:", paste(dim(expr_12L12D), collapse = " x "), "\n")

## ---------- 12L12D → LL ----------

ldll_res <- load_condition_PPRclock("data_clean/MBQN_12L12D_to_LL_PPRclock.csv")
annot_12L12D_LL <- ldll_res$annot
expr_12L12D_LL  <- ldll_res$expr

cat("[12L12D→LL] expr_12L12D_LL:", paste(dim(expr_12L12D_LL), collapse = " x "), "\n")

## ---------- 16L8D (Long Day) ----------

ldlong_res <- load_condition_PPRclock("data_clean/MBQN_16L8D_PPRclock.csv")
annot_16L8D <- ldlong_res$annot
expr_16L8D  <- ldlong_res$expr

cat("[16L8D] expr_16L8D:", paste(dim(expr_16L8D), collapse = " x "), "\n")

## ---------- 8L16D (Short Day) ----------

sd_res <- load_condition_PPRclock("data_clean/MBQN_8L16D_PPRclock.csv")
annot_8L16D <- sd_res$annot
expr_8L16D  <- sd_res$expr

cat("[8L16D] expr_8L16D:", paste(dim(expr_8L16D), collapse = " x "), "\n")

## ---------- Optional quick sanity check with a known gene ----------

test_gene <- "AT5G61380"  # TOC1

cat("Is", test_gene, "in 12L12D?      ", test_gene %in% rownames(expr_12L12D), "\n")
cat("Is", test_gene, "in 12L12D→LL?  ", test_gene %in% rownames(expr_12L12D_LL), "\n")
cat("Is", test_gene, "in 16L8D?      ", test_gene %in% rownames(expr_16L8D), "\n")
cat("Is", test_gene, "in 8L16D?      ", test_gene %in% rownames(expr_8L16D), "\n")