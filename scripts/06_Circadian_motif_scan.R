## ============================================================
## 06_Circadian_motif_scan.R
##
## Purpose
##   (1) Scan RSAT promoter FASTA files (one per photoperiod) for a small
##       panel of short circadian motifs using exact string matching on
##       both strands.
##   (2) Write per-condition motif-count tables.
##   (3) Build a single non-redundant AGI table using MAX counts across
##       conditions (promoter sequence is identical, so repeated scans are
##       redundant; MAX is a robust guard against RSAT header/entry quirks).
##   (4) For TRUE PPR oscillators (ROBUST + CANDIDATE), plot:
##       A) # motif types present (0–6)
##       B) total motif hits per gene
##
## Inputs
##   - data_raw/rsat/<COND>_promoters.txt
##   - tables/final_consensus_table.csv
##   - data_clean/PPR_only_genelist.txt
##
## Outputs
##   - tables/motif_scan_per_photoperiod/motif_counts_<COND>.csv
##   - tables/motif_scan_per_photoperiod/motif_counts_ALL_uniqueAGI_MAX.csv
##   - figures/Fig_motif_burden_TRUEoscillators_side_by_side.png
##   - tables/motif_burden_TRUEoscillators_optionA_types.csv
##   - tables/motif_burden_TRUEoscillators_optionB_totalhits.csv
##
## Notes
##   - Exact matches only: Biostrings::vcountPattern(..., fixed=TRUE)
##   - Both strands: forward + reverse complement; palindromes counted once
##   - AGI parsed from RSAT headers with (AT[1-5]G\\d{5})
## ============================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

## ----------------------------
## Settings
## ----------------------------
promoter_dir <- "data_raw/rsat"
out_dir      <- "tables/motif_scan_per_photoperiod"
fig_dir      <- "figures"

conditions <- c("12L12D", "12L12D_LL", "16L8D", "8L16D")
promoter_files <- paste0(conditions, "_promoters.txt")

motifs <- c(
  AAATATCT = "AAATATCT",
  AATATCT  = "AATATCT",
  AAAAATCT = "AAAAATCT",
  AACCAC   = "AACCAC",
  CCACAC   = "CCACAC",
  CACGTG   = "CACGTG"
)
motif_cols <- names(motifs)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## ----------------------------
## Helpers
## ----------------------------

extract_agi <- function(headers) {
  m <- stringr::str_match(headers, "(AT[1-5]G\\d{5})")
  toupper(m[, 2])
}

count_motif_both_strands <- function(seqs, motif) {
  pat <- DNAString(toupper(motif))
  rc  <- reverseComplement(pat)
  
  c_fwd <- Biostrings::vcountPattern(pat, seqs, fixed = TRUE)
  
  # Palindrome guard
  if (as.character(rc) == as.character(pat)) return(as.integer(c_fwd))
  
  c_rev <- Biostrings::vcountPattern(rc, seqs, fixed = TRUE)
  as.integer(c_fwd + c_rev)
}

scan_one_condition_counts <- function(path_fasta, cond, motifs) {
  seqs    <- Biostrings::readDNAStringSet(path_fasta)
  headers <- names(seqs)
  agi     <- extract_agi(headers)
  
  if (any(is.na(agi))) {
    warning(
      "Some FASTA headers missing AGI in ", basename(path_fasta),
      ". Examples:\n", paste(head(headers[is.na(agi)], 3), collapse = "\n")
    )
  }
  
  counts_mat <- sapply(motifs, function(m) count_motif_both_strands(seqs, m))
  counts_df  <- as.data.frame(counts_mat, stringsAsFactors = FALSE)
  
  counts_df$AGI <- agi
  counts_df$total_hits_all <- rowSums(counts_df[, names(motifs), drop = FALSE], na.rm = TRUE)
  
  gene_counts <- counts_df %>%
    filter(!is.na(AGI)) %>%
    group_by(AGI) %>%
    summarise(
      n_entries = n(),
      across(all_of(names(motifs)), ~ sum(.x, na.rm = TRUE)),
      total_hits_all = sum(total_hits_all, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(total_hits_all), AGI)
  
  out_counts <- file.path(out_dir, paste0("motif_counts_", cond, ".csv"))
  write.csv(gene_counts, out_counts, row.names = FALSE)
  
  message(sprintf("Wrote %-35s (genes=%d)", basename(out_counts), nrow(gene_counts)))
  invisible(gene_counts)
}

## ============================================================
## 1) Scan promoters per condition
## ============================================================
written <- character(0)

for (i in seq_along(conditions)) {
  cond <- conditions[i]
  fn   <- promoter_files[i]
  path <- file.path(promoter_dir, fn)
  stopifnot(file.exists(path))
  
  message("\nScanning: ", cond, " | file: ", fn)
  scan_one_condition_counts(path, cond, motifs)
  written <- c(written, paste0("motif_counts_", cond, ".csv"))
}

message("\nDone. Per-condition motif tables written to: ", out_dir)
print(written)

## ============================================================
## 2) Build ONE non-redundant table across conditions (MAX per AGI)
## ============================================================
motif_files <- file.path(out_dir, paste0("motif_counts_", conditions, ".csv"))

all_counts <- lapply(motif_files, function(f) {
  stopifnot(file.exists(f))
  x <- read.csv(f, stringsAsFactors = FALSE)
  x$AGI <- toupper(trimws(x$AGI))
  x$source_condition <- gsub("^motif_counts_|\\.csv$", "", basename(f))
  x
}) %>%
  bind_rows()

all_unique <- all_counts %>%
  group_by(AGI) %>%
  summarise(
    n_conditions_present = n_distinct(source_condition),
    n_entries_max = suppressWarnings(max(n_entries, na.rm = TRUE)),
    across(all_of(motif_cols), ~ suppressWarnings(max(.x, na.rm = TRUE))),
    .groups = "drop"
  ) %>%
  mutate(total_hits_all = rowSums(across(all_of(motif_cols)), na.rm = TRUE)) %>%
  arrange(desc(total_hits_all), AGI)

out_all <- file.path(out_dir, "motif_counts_ALL_uniqueAGI_MAX.csv")
write.csv(all_unique, out_all, row.names = FALSE)

message("\nWrote combined unique-AGI motif table (MAX across conditions):")
message("  ", out_all)
message("  Unique AGIs: ", nrow(all_unique))

## ============================================================
## 3) Plot motif burden for TRUE PPR oscillators only
## ============================================================
final_file <- "tables/final_consensus_table.csv"
ppr_file   <- "data_clean/PPR_only_genelist.txt"
stopifnot(file.exists(final_file), file.exists(ppr_file), file.exists(out_all))

ppr_only <- readLines(ppr_file, warn = FALSE) |>
  toupper() |>
  trimws()
ppr_only <- ppr_only[nzchar(ppr_only)]

true_agi <- read_csv(final_file, show_col_types = FALSE) %>%
  mutate(AGI = toupper(trimws(AGI))) %>%
  filter(consensus_candidate == TRUE) %>%     # TRUE = ROBUST + CANDIDATE in your pipeline
  filter(AGI %in% ppr_only) %>%
  distinct(AGI) %>%
  pull(AGI)

mot <- read_csv(out_all, show_col_types = FALSE) %>%
  mutate(AGI = toupper(trimws(AGI))) %>%
  filter(AGI %in% true_agi)

## ---- Option A: # motif types present (0–6) ----
tab_A <- mot %>%
  mutate(n_motif_types = rowSums(across(all_of(motif_cols), ~ as.integer(.x > 0)), na.rm = TRUE)) %>%
  count(n_motif_types, name = "n_genes") %>%
  arrange(n_motif_types)

pA <- ggplot(tab_A, aes(x = n_motif_types, y = n_genes)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, 6, by = 2), limits = c(-0.5, 6.5)) +
  scale_y_continuous(breaks = function(lims) seq(0, ceiling(lims[2] / 5) * 5, by = 5)) +
  labs(
    title = "A. Motif-type burden",
    x = "Motif types present (0–6)",
    y = "Number of TRUE oscillators"
  ) +
  theme_bw(base_size = 12)

## ---- Option B: total motif hits per gene ----
tab_B <- mot %>%
  mutate(total_hits = rowSums(across(all_of(motif_cols), ~ as.integer(.x)), na.rm = TRUE)) %>%
  count(total_hits, name = "n_genes") %>%
  arrange(total_hits)

max_hits <- max(tab_B$total_hits, na.rm = TRUE)

pB <- ggplot(tab_B, aes(x = total_hits, y = n_genes)) +
  geom_col() +
  scale_x_continuous(breaks = seq(0, max_hits, by = 2), limits = c(-0.5, max_hits + 0.5)) +
  scale_y_continuous(breaks = function(lims) seq(0, ceiling(lims[2] / 5) * 5, by = 5)) +
  labs(
    title = "B. Total motif hits",
    x = "Total motif hits per gene",
    y = "Number of TRUE oscillators"
  ) +
  theme_bw(base_size = 12)

p_combined <- pA + pB + plot_layout(ncol = 2, widths = c(1, 1))

fig_out <- file.path(fig_dir, "Fig_motif_burden_TRUEoscillators_side_by_side.png")
ggsave(fig_out, p_combined, width = 12, height = 4.8, dpi = 300)

write_csv(tab_A, "tables/motif_burden_TRUEoscillators_optionA_types.csv")
write_csv(tab_B, "tables/motif_burden_TRUEoscillators_optionB_totalhits.csv")

message("\nSaved figure:\n  ", fig_out)

print(p_combined)
