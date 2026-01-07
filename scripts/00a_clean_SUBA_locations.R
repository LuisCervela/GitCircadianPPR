## ============================================================
## 00a_clean_SUBA_locations.R
##
## Input:
##   data_raw/suba_locations.xlsx   (full SUBA export)
##
## Output:
##   data_clean/suba_location_consensus_clean.csv
##   data_clean/suba_location_consensus_clean.xlsx
##
## What it does:
##   - reads the SUBA Excel export (first sheet by default)
##   - extracts AGI locus (ATxGxxxxx) and "location consensus"
##   - standardizes AGIs, trims text, drops empty rows
##   - collapses duplicates per AGI (keeps the most informative entry)
## ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(writexl)
})

infile  <- "data_raw/suba_locations.xlsx"
out_csv <- "data_clean/suba_location_consensus_clean.csv"
out_xlsx <- "data_clean/suba_location_consensus_clean.xlsx"

stopifnot(file.exists(infile))
dir.create("data_clean", showWarnings = FALSE, recursive = TRUE)

.clean_agi <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- sub("\\..*$", "", x)                    # drop transcript suffixes if any
  m <- str_match(x, "(AT[1-5]G\\d{5})")[,2]    # extract canonical AGI pattern
  m
}

.pick_best_location <- function(loc_vec) {
  # Prefer longer, more informative strings; deprioritize blanks/NA.
  loc_vec <- unique(trimws(loc_vec))
  loc_vec <- loc_vec[nzchar(loc_vec)]
  if (length(loc_vec) == 0) return(NA_character_)
  loc_vec[which.max(nchar(loc_vec))]
}

# ---- Read Excel (sheet 1 by default) ----
raw <- readxl::read_xlsx(infile, sheet = 1, col_names = TRUE)

if (nrow(raw) == 0) stop("No rows read from: ", infile)

# ---- Identify AGI column robustly ----
cn <- colnames(raw)

agi_col <- NA_character_
# common names first
cands_agi <- c("AGI", "Locus", "Gene", "Gene ID", "GeneID", "locus", "locus_id")
agi_col <- intersect(cands_agi, cn)[1]

# fallback: search any column containing ATxGxxxxx
if (is.na(agi_col)) {
  hit_counts <- sapply(raw, function(v) sum(!is.na(.clean_agi(v))))
  if (max(hit_counts) > 0) agi_col <- names(hit_counts)[which.max(hit_counts)]
}

if (is.na(agi_col)) {
  stop("Could not infer the AGI column. Please rename the locus column to 'AGI' or 'Locus'.")
}

# ---- Identify location consensus column ----
# You said it's column 10; we still try by name first to be safe.
loc_col <- NA_character_
name_hits <- grep("consensus.*loc|loc.*consensus|location.*consensus|SUBA.*consensus",
                  cn, ignore.case = TRUE, value = TRUE)

if (length(name_hits) > 0) {
  loc_col <- name_hits[1]
} else {
  if (ncol(raw) < 10) {
    stop("Input has <10 columns, but you expected location consensus in column 10.")
  }
  loc_col <- cn[10]
}

# ---- Build clean table ----
clean <- raw %>%
  transmute(
    AGI = .clean_agi(.data[[agi_col]]),
    location_consensus = trimws(as.character(.data[[loc_col]]))
  ) %>%
  filter(!is.na(AGI), nzchar(AGI)) %>%
  mutate(
    location_consensus = na_if(location_consensus, ""),
    location_consensus = ifelse(is.na(location_consensus), NA_character_,
                                str_squish(location_consensus))
  )

if (nrow(clean) == 0) {
  stop("After cleaning, zero AGIs remained. Check that the AGI column contains ATxGxxxxx IDs.")
}

# ---- Collapse duplicates per AGI (SUBA exports sometimes repeat genes) ----
clean_dedup <- clean %>%
  group_by(AGI) %>%
  summarise(
    location_consensus = .pick_best_location(location_consensus),
    n_rows_in_input = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(AGI)

# ---- Quick diagnostics ----
cat("\nSUBA clean export diagnostics:\n")
cat("  Input rows:          ", nrow(raw), "\n")
cat("  Clean AGI rows:      ", nrow(clean), "\n")
cat("  Unique AGIs:         ", nrow(clean_dedup), "\n")
cat("  Missing locations:   ", sum(is.na(clean_dedup$location_consensus)), "\n\n")

# ---- Write outputs ----
write.csv(clean_dedup, out_csv, row.names = FALSE)
writexl::write_xlsx(clean_dedup, out_xlsx)

cat("Wrote:\n  - ", out_csv, "\n  - ", out_xlsx, "\n", sep = "")
