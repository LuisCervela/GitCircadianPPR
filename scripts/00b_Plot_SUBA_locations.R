## ============================================================
## 00b_Plot_SUBA_locations.R
## ============================================================


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ComplexUpset)
  library(ggplot2)
})

# ---------------------------
# 1) Load cleaned SUBA table
# ---------------------------
suba <- read_csv("data_clean/suba_location_consensus_clean.csv", show_col_types = FALSE) %>%
  mutate(
    AGI = toupper(trimws(AGI)),
    location_consensus = tolower(trimws(location_consensus)),
    location_consensus = str_replace_all(location_consensus, "\\s*;\\s*", ","),
    location_consensus = str_replace_all(location_consensus, "\\s*,\\s*", ",")
  )

# ---------------------------
# 2) Expand to gene × location
# ---------------------------
suba_long <- suba %>%
  separate_rows(location_consensus, sep = ",") %>%
  mutate(loc = str_trim(location_consensus)) %>%
  filter(!is.na(loc), loc != "") %>%
  distinct(AGI, loc)

# ---------------------------
# 3) Map to bins
# ---------------------------
map_loc_bin <- function(x) {
  x <- str_trim(tolower(x))
  case_when(
    x %in% c("mitochondrion", "mitochondria") ~ "mitochondrion",
    x %in% c("plastid", "chloroplast")       ~ "chloroplast",
    x %in% c("nucleus", "nuclear")           ~ "nucleus",
    x %in% c("cytosol", "cytoplasm")         ~ "cytosol",
    TRUE                                     ~ "other"
  )
}

suba_bins <- suba_long %>%
  mutate(bin = map_loc_bin(loc)) %>%
  distinct(AGI, bin)

# ---------------------------
# 4) Wide binary membership matrix (TRUE/FALSE)
# ---------------------------
bins_to_use <- c("mitochondrion", "chloroplast", "nucleus", "cytosol", "other")

suba_bin <- suba_bins %>%
  mutate(present = TRUE) %>%
  pivot_wider(
    names_from = bin,
    values_from = present,
    values_fill = FALSE
  ) %>%
  # ensure all expected bins exist as columns
  mutate(across(all_of(setdiff(bins_to_use, names(.))), ~ FALSE)) %>%
  select(AGI, all_of(bins_to_use)) %>%
  mutate(n_bins = rowSums(across(all_of(bins_to_use))))

# ---------------------------
# 5) Summaries
# ---------------------------

# genes per bin
counts_per_bin <- suba_bin %>%
  summarise(across(all_of(bins_to_use), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "bin", values_to = "n_genes") %>%
  arrange(desc(n_genes))

print(counts_per_bin)

# distribution of #bins per gene
counts_nbins <- suba_bin %>%
  count(n_bins, sort = TRUE)

print(counts_nbins)

# ---------------------------
# 6) UpSet plot (ComplexUpset syntax)
# ---------------------------
df_up <- suba_bin %>% select(AGI, all_of(bins_to_use))

p_up <- ComplexUpset::upset(
  df_up,
  intersect = bins_to_use,
  min_size = 1,
  name = "Genes"
) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  ) + 
  theme(axis.title.x = element_blank())

print(p_up)

# Optional: save
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/upset_SUBA_bins.png", p_up, width = 8, height = 5, dpi = 300)