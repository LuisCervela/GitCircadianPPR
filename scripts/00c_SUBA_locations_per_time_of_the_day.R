## ============================================================
## 00c_SUBA_locations_per_time_of_the_day.R
## ============================================================
library(dplyr)
library(readr)
library(tidyr)
library(writexl)

## ----------------------------
## 1) Phase binning (4h windows)
## ----------------------------
phase_to_zt_hour <- function(x) {
  # round to nearest hour, wrap to 0..23
  as.integer(round(x) %% 24)
}

zt_to_bin4h <- function(zt_hour) {
  # your bins:
  # Late night: 0–3
  # Morning: 4–7
  # Midday: 8–11
  # Afternoon: 12–15
  # Evening: 16–19
  # Late evening: 20–23
  cut(
    zt_hour,
    breaks = c(-1, 3, 7, 11, 15, 19, 23),
    labels = c("Dawn/early day (ZT0–3)",
               "Morning (ZT4–7)",
               "Late day (ZT8–11)",
               "Early night (ZT12–15)",
               "Mid night (ZT16–19)",
               "Late night/pre-dawn (ZT20–23)"),
    right = TRUE
  )
}

bin_levels <- c("Dawn/early day (ZT0–3)",
                "Morning (ZT4–7)",
                "Late day (ZT8–11)",
                "Early night (ZT12–15)",
                "Mid night (ZT16–19)",
                "Late night/pre-dawn (ZT20–23)")

## ----------------------------
## 2) Load TRUE oscillators (robust + candidate) with phases
## ----------------------------
final <- read_csv("tables/final_consensus_table.csv", show_col_types = FALSE) %>%
  mutate(
    AGI = toupper(trimws(AGI)),
    condition = as.character(condition),
    phase_est = suppressWarnings(as.numeric(phase_est))
  ) %>%
  filter(consensus_candidate == TRUE, is.finite(phase_est)) %>%  # TRUE = robust + candidate (your pipeline)
  mutate(
    zt_hour = phase_to_zt_hour(phase_est),
    phase_bin4h = factor(zt_to_bin4h(zt_hour), levels = bin_levels)
  )

ppr_only <- read_lines("data_clean/PPR_only_genelist.txt") %>%
  toupper() %>% trimws() %>% (\(x) x[nzchar(x)])()

final <- final %>%
  filter(AGI %in% ppr_only)

## ----------------------------
## 3) Load SUBA memberships
## You said you already have a cleaned file like:
##   data_clean/suba_location_consensus_clean.csv
## If instead you have a different path, change it here.
## ----------------------------
suba <- read_csv("data_clean/suba_location_consensus_clean.csv", show_col_types = FALSE) %>%
  mutate(
    AGI = toupper(trimws(AGI)),
    location_consensus = tolower(trimws(location_consensus))
  )

# Expand multi-localizations -> one row per AGI x raw location token
suba_long <- suba %>%
  mutate(location_consensus = gsub(";", ",", location_consensus),
         location_consensus = gsub("\\s*,\\s*", ",", location_consensus)) %>%
  separate_rows(location_consensus, sep = ",") %>%
  mutate(loc = trimws(location_consensus)) %>%
  filter(!is.na(loc), loc != "") %>%
  distinct(AGI, loc)

# Map raw SUBA strings into the 5 bins you were using internally
map_loc5 <- function(x) {
  x <- trimws(tolower(x))
  dplyr::case_when(
    x %in% c("mitochondrion", "mitochondria") ~ "mito",
    x %in% c("plastid", "chloroplast")       ~ "chloro",
    x %in% c("nucleus", "nuclear")           ~ "nucleus",
    x %in% c("cytosol", "cytoplasm")         ~ "cytosol",
    TRUE                                     ~ "other"
  )
}

suba_sets <- suba_long %>%
  mutate(bin = map_loc5(loc)) %>%
  distinct(AGI, bin)

# Binary membership (one row per gene)
suba_bin <- suba_sets %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = bin, values_from = present, values_fill = 0L) %>%
  mutate(
    mito    = as.integer(mito > 0),
    chloro  = as.integer(chloro > 0),
    nucleus = as.integer(nucleus > 0),
    cytosol = as.integer(cytosol > 0),
    other   = as.integer(other > 0)
  ) %>%
  select(AGI, mito, chloro, nucleus, cytosol, other)

## ----------------------------
## 4) Strict 4-class location assignment
## dual = mito+chloro ONLY (no nucleus/cytosol/other)
## ----------------------------
suba_loc4 <- suba_bin %>%
  mutate(
    loc4 = case_when(
      mito == 1 & chloro == 1 & nucleus == 0 & cytosol == 0 & other == 0 ~ "dual",
      mito == 1 & chloro == 0                                           ~ "mito",
      chloro == 1 & mito == 0                                           ~ "chloro",
      TRUE                                                              ~ "other"
    )
  ) %>%
  select(AGI, loc4)

## ----------------------------
## 5) Build ONE wide table per condition
## rows = phase bins, cols = mito/chloro/dual/other
## ----------------------------
make_phase_location_table <- function(cond_name) {
  df <- final %>%
    filter(condition == cond_name) %>%
    left_join(suba_loc4, by = "AGI") %>%
    mutate(loc4 = ifelse(is.na(loc4), "other", loc4)) %>%
    count(phase_bin4h, loc4, name = "n_genes") %>%
    tidyr::complete(
      phase_bin4h = factor(bin_levels, levels = bin_levels),
      loc4 = c("mito","chloro","dual","other"),
      fill = list(n_genes = 0L)
    ) %>%
    pivot_wider(names_from = loc4, values_from = n_genes) %>%
    arrange(phase_bin4h)
  
  df
}

conds <- c("12L12D","12L12D_LL","16L8D","8L16D")
tables_by_condition <- setNames(lapply(conds, make_phase_location_table), conds)

# Write each as CSV (optional)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
for (nm in names(tables_by_condition)) {
  write_csv(tables_by_condition[[nm]], file.path("tables", paste0("TRUE_PPR_phasebin4h_loc4_", nm, ".csv")))
}

# Write all into one Excel file (one sheet per condition)
write_xlsx(tables_by_condition, path = "tables/TRUE_PPR_phasebin4h_loc4_by_condition.xlsx")


tables_by_condition[["12L12D_LL"]] # view example
