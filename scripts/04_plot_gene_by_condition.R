library(dplyr)
library(ggplot2)

plot_gene_condition <- function(gene_id,
                                expr_mat,
                                annot,
                                condition_label = "condition") {
  if (!"column_name" %in% colnames(annot)) {
    stop("Annotation data.frame must contain a 'column_name' column.")
  }
  if (!all(annot$column_name %in% colnames(expr_mat))) {
    stop("Some annot$column_name are not columns in expr_mat.")
  }
  if (!gene_id %in% rownames(expr_mat)) {
    stop(paste("Gene", gene_id, "not found in expression matrix."))
  }
  if (!all(c("gear_id", "time") %in% colnames(annot))) {
    stop("Annotation must contain 'gear_id' and 'time' columns.")
  }
  
  gene_values <- expr_mat[gene_id, ]
  
  gene_df <- data.frame(
    AGI         = gene_id,
    expression  = as.numeric(gene_values),
    column_name = names(gene_values),
    stringsAsFactors = FALSE
  ) %>%
    left_join(annot, by = "column_name")
  
  gene_summary <- gene_df %>%
    group_by(gear_id, time) %>%
    summarise(
      expression_mean = mean(expression, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  p_all <- ggplot(gene_summary,
                  aes(x = time,
                      y = expression_mean,
                      color = factor(gear_id),
                      group = gear_id)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = paste0(gene_id, " • ", condition_label, " across datasets"),
      x = "Time",
      y = "MBQN expression (mean per time)",
      color = "GEAR ID"
    ) +
    theme_bw(base_size = 14)
  
  p_facets <- ggplot(gene_summary,
                     aes(x = time,
                         y = expression_mean,
                         group = 1)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(size = 2, color = "black") +
    facet_wrap(~ gear_id, scales = "free_x") +
    labs(
      title = paste0(gene_id, " • ", condition_label, " per dataset"),
      x = "Time",
      y = "MBQN expression (mean per time)"
    ) +
    theme_bw(base_size = 14)
  
  list(
    data        = gene_summary,
    plot_all    = p_all,
    plot_panels = p_facets
  )
}


get_condition_data <- function(condition) {
  # condition is one of: "12L12D", "12L12D_LL", "16L8D", "8L16D"
  cond <- match.arg(condition,
                    choices = c("12L12D", "12L12D_LL", "16L8D", "8L16D"))
  
  if (cond == "12L12D") {
    if (!exists("expr_12L12D") || !exists("annot_12L12D")) {
      stop("expr_12L12D / annot_12L12D not found. Did you run 03_load_PPRclock_all_conditions.R?")
    }
    return(list(expr = expr_12L12D,
                annot = annot_12L12D,
                label = "LD (12L12D)"))
  }
  
  if (cond == "12L12D_LL") {
    if (!exists("expr_12L12D_LL") || !exists("annot_12L12D_LL")) {
      stop("expr_12L12D_LL / annot_12L12D_LL not found.")
    }
    return(list(expr = expr_12L12D_LL,
                annot = annot_12L12D_LL,
                label = "LD→LL (12L12D→LL)"))
  }
  
  if (cond == "16L8D") {
    if (!exists("expr_16L8D") || !exists("annot_16L8D")) {
      stop("expr_16L8D / annot_16L8D not found.")
    }
    return(list(expr = expr_16L8D,
                annot = annot_16L8D,
                label = "Long Day (16L8D)"))
  }
  
  if (cond == "8L16D") {
    if (!exists("expr_8L16D") || !exists("annot_8L16D")) {
      stop("expr_8L16D / annot_8L16D not found.")
    }
    return(list(expr = expr_8L16D,
                annot = annot_8L16D,
                label = "Short Day (8L16D)"))
  }
}

plot_gene_by_condition <- function(gene_id,
                                   condition = c("12L12D", "12L12D_LL", "16L8D", "8L16D")) {
  cond <- match.arg(condition)
  dat  <- get_condition_data(cond)
  
  plot_gene_condition(
    gene_id        = gene_id,
    expr_mat       = dat$expr,
    annot          = dat$annot,
    condition_label = dat$label
  )
}


gene_id   <- "AT2G46830"      # CCA1
condition <- "12L12D"         # options: "12L12D", "12L12D_LL", "16L8D", "8L16D"

res <- plot_gene_by_condition(gene_id, condition)

# All datasets overlayed
res$plot_all

# One panel per GEAR ID
res$plot_panels
