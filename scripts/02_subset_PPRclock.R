## 02_subset_PPRclock.R
## Subset MBQN condition tables to PPR + clock gene list,
## and create both full (with metadata) and expression-only versions.

library(data.table)

## ---------- 1. Read gene list ----------

gene_list_file <- "data_clean/PPRClockgenelist.md"

genes_df <- read.table(
  gene_list_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

genes_vec <- unique(genes_df[[1]])

# drop possible header variants
genes_vec <- genes_vec[!(genes_vec %in% c("AGI", "agi"))]

## ---------- 2. Helper: subset one file by gene list ----------

subset_by_genes <- function(infile, outfile, genes_vec) {
  dat <- fread(infile, data.table = FALSE, check.names = FALSE)
  
  # First 3 rows = metadata (Time, Accession, Sample)
  meta <- dat[1:3, ]
  
  # Expression rows
  expr <- dat[-c(1:3), ]
  colnames(expr)[1] <- "AGI"
  
  # Clean gene IDs
  genes_vec_clean <- toupper(trimws(gsub("^[-*]\\s*", "", genes_vec)))
  expr$AGI_clean  <- toupper(trimws(sub("\\..*$", "", expr$AGI)))
  
  # Subset using cleaned IDs
  expr_sub <- expr[expr$AGI_clean %in% genes_vec_clean, ]
  
  # Drop helper and restore original first column name
  expr_sub$AGI_clean <- NULL
  colnames(expr_sub)[1] <- colnames(dat)[1]
  
  # Rebuild table: metadata + filtered expression
  dat_sub <- rbind(meta, expr_sub)
  
  fwrite(dat_sub, outfile)
  cat("Wrote", nrow(expr_sub), "genes to", outfile, "\n")
}

## ---------- 3. Helper: make expression-only version ----------

make_expr_only <- function(infile, outfile) {
  x <- fread(infile, data.table = FALSE, check.names = FALSE)
  x_expr <- x[-c(1:3), ]
  colnames(x_expr)[1] <- "AGI"
  fwrite(x_expr, outfile)
  cat("Saved expression-only:", outfile, "\n")
}

## ---------- 4. Apply to each condition ----------

## 12L12D (LD)
subset_by_genes(
  infile  = "data_clean/MBQN_12L12D.csv",
  outfile = "data_clean/MBQN_12L12D_PPRclock.csv",
  genes_vec = genes_vec
)
make_expr_only(
  infile  = "data_clean/MBQN_12L12D_PPRclock.csv",
  outfile = "data_clean/MBQN_12L12D_PPRclock_expr_only.csv"
)

## 12L12D → LL
subset_by_genes(
  infile  = "data_clean/MBQN_12L12D_to_LL.csv",
  outfile = "data_clean/MBQN_12L12D_to_LL_PPRclock.csv",
  genes_vec = genes_vec
)
make_expr_only(
  infile  = "data_clean/MBQN_12L12D_to_LL_PPRclock.csv",
  outfile = "data_clean/MBQN_12L12D_to_LL_PPRclock_expr_only.csv"
)

## 16L8D (Long Day)
subset_by_genes(
  infile  = "data_clean/MBQN_16L8D.csv",
  outfile = "data_clean/MBQN_16L8D_PPRclock.csv",
  genes_vec = genes_vec
)
make_expr_only(
  infile  = "data_clean/MBQN_16L8D_PPRclock.csv",
  outfile = "data_clean/MBQN_16L8D_PPRclock_expr_only.csv"
)

## 8L16D (Short Day)
subset_by_genes(
  infile  = "data_clean/MBQN_8L16D.csv",
  outfile = "data_clean/MBQN_8L16D_PPRclock.csv",
  genes_vec = genes_vec
)
make_expr_only(
  infile  = "data_clean/MBQN_8L16D_PPRclock.csv",
  outfile = "data_clean/MBQN_8L16D_PPRclock_expr_only.csv"
)


