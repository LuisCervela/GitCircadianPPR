library(data.table)

## ---------- 1. Read raw MBQN file ----------

raw_file <- "data_raw/MBQN_20250715.csv"
if (!file.exists(raw_file)) {
  stop(
    "Missing: ", raw_file, "\n",
    "Download MBQN from GEAR (https://scientist.xsrv.jp/wp-content/uploads/2025/07/MBQN_20250715.csv) ",
    "or from Zenodo (DOI: 10.5281/zenodo.15910335), ",
    "then place it in data_raw/ with the exact filename."
  )
}



mbqn <- fread(
  "data_raw/MBQN_20250715.csv",
  header      = TRUE,
  sep         = ",",
  data.table  = FALSE,
  check.names = FALSE
)

## ---------- 2. Select only GEAR IDs of interest ----------

cn <- colnames(mbqn)

gear_id_from_col <- suppressWarnings(as.numeric(sub("\\..*$", "", cn)))

keep_ids <- c(
  4, 5, 8, 14, 15, 17,
  26, 27, 28, 29, 30, 31, 32,
  33, 34, 35, 46, 48, 50, 51, 52,
  58, 75, 76, 77, 78
)

cols_keep <- which(
  gear_id_from_col %in% keep_ids |
    cn == "GEAR ID"
)

mbqn_trimmed <- mbqn[, cols_keep]

fwrite(
  mbqn_trimmed,
  file = "data_clean/MBQN_trimmed_selected_GEARS.csv"
)

## ---------- 3. Build sample annotation ----------

cn_trim <- colnames(mbqn_trimmed)
gear_id_trim <- suppressWarnings(as.numeric(sub("\\..*$", "", cn_trim)))

time_vec      <- as.character(mbqn_trimmed[1, ])
accession_vec <- as.character(mbqn_trimmed[2, ])
sample_vec    <- as.character(mbqn_trimmed[3, ])

sample_annotation <- data.frame(
  column_name = cn_trim,
  gear_id     = gear_id_trim,
  time        = time_vec,
  accession   = accession_vec,
  sample_name = sample_vec,
  stringsAsFactors = FALSE
)

write.csv(
  sample_annotation,
  "data_clean/sample_annotation.csv",
  row.names = FALSE
)

## ---------- 4. Define photoperiod groups ----------

ids_12L12D    <- c(8, 14, 17, 58, 78)
ids_12L12D_LL <- c(4, 5, 26, 27, 75, 76, 77)
ids_16L8D     <- c(28, 29, 30, 31, 32, 51)
ids_8L16D     <- c(15, 33, 34, 35, 46, 48, 50, 52)

gear_col_idx <- which(cn_trim == "GEAR ID")

cols_12L12D    <- c(gear_col_idx, which(gear_id_trim %in% ids_12L12D))
cols_12L12D_LL <- c(gear_col_idx, which(gear_id_trim %in% ids_12L12D_LL))
cols_16L8D     <- c(gear_col_idx, which(gear_id_trim %in% ids_16L8D))
cols_8L16D     <- c(gear_col_idx, which(gear_id_trim %in% ids_8L16D))

## ---------- 5. Split trimmed MBQN into condition-specific matrices ----------

mbqn_12L12D    <- mbqn_trimmed[, cols_12L12D]
mbqn_12L12D_LL <- mbqn_trimmed[, cols_12L12D_LL]
mbqn_16L8D     <- mbqn_trimmed[, cols_16L8D]
mbqn_8L16D     <- mbqn_trimmed[, cols_8L16D]

fwrite(mbqn_12L12D,    "data_clean/MBQN_12L12D.csv")
fwrite(mbqn_12L12D_LL, "data_clean/MBQN_12L12D_LL.csv")
fwrite(mbqn_16L8D,     "data_clean/MBQN_16L8D.csv")
fwrite(mbqn_8L16D,     "data_clean/MBQN_8L16D.csv")
