## ============================================================
## 02_subset_PPRclock.R
##
## Purpose
##   - Read a gene list (PPR + clock genes) from a simple AGI list file.
##   - Subset each photoperiod MBQN table (from Script 01) to that gene list.
##   - Write, for each condition:
##       (1) a "full" table (metadata rows + expression rows)
##       (2) an "expression-only" table (expression rows only)
##   - Additionally, split the master list into:
##       - clock-only gene list
##       - PPR-only gene list  (= all genes minus clock genes)
##
## Inputs (from Script 01)
##   - data_clean/MBQN_12L12D.csv
##   - data_clean/MBQN_12L12D_LL.csv
##   - data_clean/MBQN_16L8D.csv
##   - data_clean/MBQN_8L16D.csv
##
## Input gene list
##   - data_clean/PPRClockgenelist.md   (AGI codes, one per line; markdown bullets OK)
##
## Outputs (written to data_clean/)
##   - MBQN_<COND>_PPRclock.csv
##   - MBQN_<COND>_PPRclock_expr_only.csv
##   - PPR_only_genelist.txt
##   - clock_only_genelist.txt
##
## Notes
##   - MBQN condition tables are assumed to have 3 metadata rows:
##       row 1 = time, row 2 = accession, row 3 = sample name
##   - Gene IDs are matched case-insensitively and ignoring isoform suffix (e.g., ".1")
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

## ----------------------------
## 0) Paths and basic checks
## ----------------------------
indir  <- "data_clean"
outdir <- "data_clean"  # keep outputs together for now

gene_list_file <- file.path(indir, "PPRClockgenelist.md")

if (!file.exists(gene_list_file)) {
  stop("Missing gene list file: ", gene_list_file)
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ----------------------------
## 1) Read and clean AGI gene list
## ----------------------------
read_agi_list <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- toupper(trimws(x))
  x <- gsub("^[-*]\\s*", "", x)     # strip markdown bullets (e.g., "- AT1G...")
  x <- sub("\\..*$", "", x)         # drop isoform suffix if present
  x <- x[nzchar(x)]
  x <- x[grepl("^AT[1-5]G\\d{5}$", x)]  # keep valid Arabidopsis AGI codes
  unique(x)
}

all_genes <- read_agi_list(gene_list_file)

if (length(all_genes) == 0) {
  stop("No valid AGI codes were parsed from: ", gene_list_file)
}

## ----------------------------
## 2) Define clock genes and derive PPR-only list
## ----------------------------
clock_genes <- unique(c(
  "AT2G46830",
  "AT1G01060",
  "AT3G09600",
  "AT5G61380",
  "AT5G24470",
  "AT5G02810",
  "AT2G46790",
  "AT1G22770",
  "AT2G25930",
  "AT2G40080",
  "AT3G46640",
  "AT5G64170",
  "AT3G54500"
))

missing_clock <- setdiff(clock_genes, all_genes)
if (length(missing_clock) > 0) {
  warning("These clock genes were not found in ", gene_list_file, ": ",
          paste(missing_clock, collapse = ", "))
}

ppr_only <- setdiff(all_genes, clock_genes)

cat("Loaded genes from:", gene_list_file, "\n")
cat("  total genes parsed:", length(all_genes), "\n")
cat("  clock genes present:", length(intersect(clock_genes, all_genes)), "\n")
cat("  PPR-only genes:     ", length(ppr_only), "\n")

writeLines(ppr_only,    file.path(outdir, "PPR_only_genelist.txt"))
writeLines(clock_genes, file.path(outdir, "clock_only_genelist.txt"))

cat("Wrote PPR-only and clock-only gene lists to:", outdir, "\n")

## ----------------------------
## 3) Helper: subset an MBQN condition file by gene list
## ----------------------------
subset_by_genes <- function(infile, outfile, genes_vec) {
  if (!file.exists(infile)) stop("Missing input MBQN table: ", infile)
  
  dat <- fread(infile, data.table = FALSE, check.names = FALSE)
  
  if (nrow(dat) < 4) {
    stop("Input MBQN table seems too small (expected >= 4 rows): ", infile)
  }
  
  # First 3 rows = metadata
  meta <- dat[1:3, , drop = FALSE]
  
  # Expression rows
  expr <- dat[-c(1:3), , drop = FALSE]
  colnames(expr)[1] <- "AGI"
  
  # Clean gene IDs for matching
  genes_vec_clean <- toupper(trimws(gsub("^[-*]\\s*", "", genes_vec)))
  genes_vec_clean <- sub("\\..*$", "", genes_vec_clean)
  
  expr$AGI_clean <- toupper(trimws(sub("\\..*$", "", expr$AGI)))
  
  # Subset
  expr_sub <- expr[expr$AGI_clean %in% genes_vec_clean, , drop = FALSE]
  
  # Drop helper column and restore original first column name
  expr_sub$AGI_clean <- NULL
  colnames(expr_sub)[1] <- colnames(dat)[1]
  
  # Rebuild: metadata + filtered expression
  dat_sub <- rbind(meta, expr_sub)
  
  fwrite(dat_sub, outfile)
  cat("Wrote", nrow(expr_sub), "genes to", outfile, "\n")
  invisible(nrow(expr_sub))
}

## ----------------------------
## 4) Helper: make expression-only version (drop metadata)
## ----------------------------
make_expr_only <- function(infile, outfile) {
  if (!file.exists(infile)) stop("Missing input file for expr-only export: ", infile)
  
  x <- fread(infile, data.table = FALSE, check.names = FALSE)
  if (nrow(x) < 4) stop("File too small to drop metadata rows 1:3: ", infile)
  
  x_expr <- x[-c(1:3), , drop = FALSE]
  colnames(x_expr)[1] <- "AGI"
  
  fwrite(x_expr, outfile)
  cat("Saved expression-only:", outfile, "\n")
}

## ----------------------------
## 5) Apply to each condition
## ----------------------------
condition_map <- list(
  `12L12D`    = file.path(indir, "MBQN_12L12D.csv"),
  `12L12D_LL` = file.path(indir, "MBQN_12L12D_LL.csv"),
  `16L8D`     = file.path(indir, "MBQN_16L8D.csv"),
  `8L16D`     = file.path(indir, "MBQN_8L16D.csv")
)

for (cond in names(condition_map)) {
  infile <- condition_map[[cond]]
  
  out_full <- file.path(outdir, paste0("MBQN_", cond, "_PPRclock.csv"))
  out_expr <- file.path(outdir, paste0("MBQN_", cond, "_PPRclock_expr_only.csv"))
  
  subset_by_genes(infile = infile, outfile = out_full, genes_vec = all_genes)
  make_expr_only(infile = out_full, outfile = out_expr)
}

cat("Done. Subset tables written to:", outdir, "\n")


