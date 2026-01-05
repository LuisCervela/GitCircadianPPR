## ============================================================
## 07_Radar_plot_TRUE_PPR.R
##
## Purpose
##   Make 1 radar plot per photoperiod showing the distribution of peak
##   phase (phase_est) for TRUE PPR oscillators (robust + candidates).
##
##   Plot requirements implemented:
##     - ZT0 at the TOP (north) and time increases CLOCKWISE
##     - ZT labels every 2 hours, outside the circle (ZT0..ZT22; no ZT24)
##     - Internal grid lines subtle
##     - A radial axis (#genes) with tick labels
##     - A "Total TRUE PPRs: N" label at the bottom
##     - Output: PNG only
##
## Inputs
##   - tables/final_consensus_table.csv   (from Script 04 rhythmicity)
##   - data_clean/PPR_only_genelist.txt   (from Script 02)
##
## Outputs
##   - figures/radar_TRUE_PPR_peak_<COND>_polished.png
## ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

final_file <- "tables/final_consensus_table.csv"
ppr_file   <- "data_clean/PPR_only_genelist.txt"
stopifnot(file.exists(final_file), file.exists(ppr_file))

## ----------------------------
## Load inputs
## ----------------------------
ppr_only <- readLines(ppr_file, warn = FALSE) |>
  trimws() |>
  toupper()
ppr_only <- ppr_only[nzchar(ppr_only)]

final_cons <- fread(final_file, data.table = FALSE) %>%
  mutate(
    AGI = toupper(trimws(sub("\\..*$", "", AGI))),
    condition = as.character(condition),
    phase_est = suppressWarnings(as.numeric(phase_est)),
    # "TRUE oscillator" = candidate OR robust in this pipeline
    consensus_candidate = as.logical(consensus_candidate)
  )

conditions <- c("12L12D", "12L12D_LL", "16L8D", "8L16D")

## ----------------------------
## Geometry helpers
## ----------------------------
phase_to_zt_bin <- function(phase_hours) {
  ph <- phase_hours %% 24
  as.integer(round(ph)) %% 24
}

# Convert ZT (0..24) to radians with ZT0 at top and clockwise direction
zt_to_theta <- function(zt) {
  (pi/2) - (2*pi) * (zt / 24)
}

circle_df <- function(r, n = 360) {
  t <- seq(0, 2*pi, length.out = n)
  data.frame(x = r*cos(t), y = r*sin(t), r = r)
}

## ============================================================
## Main loop: one plot per condition
## ============================================================
for (cond in conditions) {
  
  df_true <- final_cons %>%
    filter(condition == cond, AGI %in% ppr_only) %>%
    filter(consensus_candidate == TRUE) %>%     # robust + candidate
    filter(is.finite(phase_est))
  
  n_total <- n_distinct(df_true$AGI)
  if (n_total == 0) {
    message("No TRUE PPRs with finite phase_est for: ", cond, " (skipping).")
    next
  }
  
  ## 1h bins
  counts <- df_true %>%
    mutate(ZT = phase_to_zt_bin(phase_est)) %>%
    count(ZT, name = "n_genes") %>%
    right_join(data.frame(ZT = 0:23), by = "ZT") %>%
    mutate(n_genes = ifelse(is.na(n_genes), 0L, n_genes)) %>%
    arrange(ZT)
  
  # close the curve by repeating ZT0 at ZT24
  counts_closed <- bind_rows(
    counts,
    counts %>% filter(ZT == 0) %>% mutate(ZT = 24)
  )
  
  max_r <- max(counts_closed$n_genes, na.rm = TRUE)
  if (!is.finite(max_r) || max_r <= 0) {
    message("max_r <= 0 for: ", cond, " (skipping).")
    next
  }
  
  df_path <- counts_closed %>%
    mutate(
      theta = zt_to_theta(ZT),
      x = n_genes * cos(theta),
      y = n_genes * sin(theta)
    )
  
  ## ---- Grid circles + spokes (subtle) ----
  r_breaks <- pretty(c(0, max_r), n = 4)
  r_breaks <- unique(r_breaks[r_breaks > 0])
  
  grid_circles <- bind_rows(lapply(r_breaks, circle_df))
  
  zt_ticks <- seq(0, 22, by = 2)  # labels: ZT0..ZT22, no ZT24
  tick_theta <- zt_to_theta(zt_ticks)
  
  spokes <- data.frame(
    x = 0, y = 0,
    xend = max_r * cos(tick_theta),
    yend = max_r * sin(tick_theta),
    ZT = zt_ticks
  )
  
  ## ZT labels outside the circle
  label_r <- max_r * 1.14
  tick_labels <- data.frame(
    x = label_r * cos(tick_theta),
    y = label_r * sin(tick_theta),
    lab = paste0("ZT", zt_ticks)
  )
  
  ## ---- Radial axis (# genes) ----
  axis_theta <- pi  # left/west
  axis_df <- data.frame(
    x = 0, y = 0,
    xend = max_r * cos(axis_theta),
    yend = max_r * sin(axis_theta)
  )
  
  radial_labels <- data.frame(
    r = r_breaks,
    x = r_breaks * cos(axis_theta) - 0.12 * max_r,
    y = r_breaks * sin(axis_theta),
    lab = as.character(r_breaks)
  )
  
  ## ---- Total label bottom (south) ----
  total_x <- 0
  total_y <- -(max_r * 1.30)
  
  p <- ggplot() +
    # subtle circles
    geom_path(
      data = grid_circles,
      aes(x = x, y = y, group = r),
      linewidth = 0.5, alpha = 0.12
    ) +
    # subtle spokes
    geom_segment(
      data = spokes,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.5, alpha = 0.12
    ) +
    # ZT labels outside
    geom_text(
      data = tick_labels,
      aes(x = x, y = y, label = lab),
      size = 4
    ) +
    # radial axis + numeric ticks
    geom_segment(
      data = axis_df,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.6, alpha = 0.30
    ) +
    geom_text(
      data = radial_labels,
      aes(x = x, y = y, label = lab),
      size = 3.6, alpha = 0.65
    ) +
    annotate(
      "text",
      x = (max_r * cos(axis_theta) - 0.22 * max_r),
      y = (max_r * sin(axis_theta) + 0.10 * max_r),
      label = "gene count",
      fontface = "bold",
      size = 4
    ) +
    # radar path
    geom_path(
      data = df_path,
      aes(x = x, y = y),
      linewidth = 1.2
    ) +
    geom_point(
      data = df_path,
      aes(x = x, y = y),
      size = 2.0
    ) +
    # total label bottom
    annotate(
      "label",
      x = total_x, y = total_y,
      label = paste0("Total PPRs: ", n_total),
      label.size = 0.25
    ) +
    coord_equal(clip = "off") +
    labs(
      title = paste0("Rythmic PPR peak-time distribution (", cond, ")"),
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 20, r = 25, b = 50, l = 35)
    )
  
  out_png <- file.path("figures", paste0("radar_TRUE_PPR_peak_", cond, "_polished.png"))
  ggsave(out_png, plot = p, width = 7, height = 7, dpi = 300)
  
  cat("Saved:", out_png, "\n")
}

cat("\nCreated:\n")
print(list.files("figures", pattern = "_polished\\.png$", full.names = FALSE))

