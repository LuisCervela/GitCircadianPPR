## ============================================================
## 06_Circadian_motif_scan.R
##
## Purpose
##   Scan RSAT-retrieved promoter FASTA files (one per photoperiod) for
##   a small set of short consensus motifs (string matches; fixed).
##
##   For each photoperiod:
##     - Count occurrences of each motif in each promoter sequence
##       on BOTH strands (forward + reverse-complement).
##     - Collapse to ONE row per AGI (sum across multiple entries if present).
##     - Write a CSV table of motif counts per gene.
##
## Inputs
##   Directory (default):
##     data_raw/rsat/
##   Files (default):
##     12L12D_promoters.txt
##     12L12D_LL_promoters.txt
##     16L8D_promoters.txt
##     8L16D_promoters.txt
##
## Outputs
##   Directory:
##     tables/motif_scan_per_photoperiod/
##   Files:
##     motif_counts_12L12D.csv
##     motif_counts_12L12D_LL.csv
##     motif_counts_16L8D.csv
##     motif_counts_8L16D.csv
##
## Notes
##   - Uses Biostrings::vcountPattern with fixed=TRUE (exact matches).
##   - Palindromic motifs (motif == reverseComplement(motif)) are counted once.
##   - AGI is parsed from RSAT FASTA headers using regex (AT[1-5]Gxxxxx).
## ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(stringr)
})

## ----------------------------
## USER INPUTS
## ----------------------------
promoter_dir <- "data_raw/rsat"

files <- c(
  "12L12D_promoters.txt",
  "12L12D_LL_promoters.txt",
  "16L8D_promoters.txt",
  "8L16D_promoters.txt"
)

motifs <- c(
  AAATATCT = "AAATATCT",
  AATATCT  = "AATATCT",
  AAAAATCT = "AAAAATCT",
  AACCAC   = "AACCAC",
  CCACAC   = "CCACAC",
  CACGTG   = "CACGTG"
)

out_dir <- "tables/motif_scan_per_photoperiod"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ----------------------------
## HELPERS
## ----------------------------

# Extract AGI codes from RSAT FASTA headers
extract_agi <- function(headers) {
  m <- stringr::str_match(headers, "(AT[1-5]G\\d{5})")
  toupper(m[, 2])
}

# Count motif occurrences on both strands (forward + reverse complement),
# without double-counting palindromes.
count_motif_both_strands <- function(seqs, motif) {
  pat <- DNAString(toupper(motif))
  rc  <- reverseComplement(pat)
  
  c_fwd <- Biostrings::vcountPattern(pat, seqs, fixed = TRUE)
  
  # Palindrome guard (e.g., CACGTG is palindromic)
  if (as.character(rc) == as.character(pat)) {
    return(as.integer(c_fwd))
  }
  
  c_rev <- Biostrings::vcountPattern(rc, seqs, fixed = TRUE)
  as.integer(c_fwd + c_rev)
}

# Scan one promoter FASTA file and write one CSV for that condition
scan_one_condition_counts <- function(path_fasta, cond, motifs, out_dir) {
  
  seqs    <- Biostrings::readDNAStringSet(path_fasta)
  headers <- names(seqs)
  agi     <- extract_agi(headers)
  
  if (any(is.na(agi))) {
    warning(
      "Some FASTA headers missing AGI in ", basename(path_fasta),
      ". Examples:\n", paste(head(headers[is.na(agi)], 3), collapse = "\n")
    )
  }
  
  # count motifs per entry
  counts_mat <- sapply(motifs, function(m) count_motif_both_strands(seqs, m))
  counts_df  <- as.data.frame(counts_mat, stringsAsFactors = FALSE)
  
  counts_df$AGI <- agi
  counts_df$total_hits_all <- rowSums(counts_df[, names(motifs), drop = FALSE], na.rm = TRUE)
  
  # collapse to one row per gene (sum across multiple entries if present)
  gene_counts <- counts_df %>%
    dplyr::filter(!is.na(AGI)) %>%
    dplyr::group_by(AGI) %>%
    dplyr::summarise(
      n_entries = dplyr::n(),
      dplyr::across(dplyr::all_of(names(motifs)), ~ sum(.x, na.rm = TRUE)),
      total_hits_all = sum(total_hits_all, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(total_hits_all), AGI)
  
  out_counts <- file.path(out_dir, paste0("motif_counts_", cond, ".csv"))
  write.csv(gene_counts, out_counts, row.names = FALSE)
  
  cat(sprintf("  Wrote %-35s  (genes=%d)\n", basename(out_counts), nrow(gene_counts)))
  invisible(gene_counts)
}

## ----------------------------
## RUN
## ----------------------------
written <- character(0)

for (fn in files) {
  path <- file.path(promoter_dir, fn)
  stopifnot(file.exists(path))
  
  # Derive condition label from filename
  cond <- sub("_promoters\\..*$", "", fn)  # e.g. 12L12D_LL
  
  cat("\nScanning:", cond, " | file:", fn, "\n")
  scan_one_condition_counts(path, cond, motifs, out_dir)
  
  written <- c(written, paste0("motif_counts_", cond, ".csv"))
}

cat("\nDone. Files written to:", out_dir, "\n")
print(written)