rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================

# ==============================================================================

library(tidyverse)
library(data.table)
library(zoo)
library(ggh4x)
library(patchwork)
library(cowplot)
library(ggridges)
library(fdapace)

# ---- Paths ----
dat_dir = paste0(getwd(),"/data/")
fig_dir = paste0(getwd(),"/figures/")
timestamp <- format(Sys.time(), "%d%m%Y")

# ==============================================================================
# 1. Read and summarize data 
# ==============================================================================


#canopy dimensions from the field
my_dat1<-read.csv(paste0(dat_dir, "04_dat_field_canopy_rvn.csv"))


# ==============================================================================
# 2. PACE analysis (platform and field)
# ==============================================================================

#function (platform)
run_fpca <- function(df, group_var, time_var, trait_col) {
  
  # Filter NAs for both trait and time variable
  df <- df %>% 
    dplyr::filter(!is.na(.data[[trait_col]]), !is.na(.data[[time_var]])) %>% 
    droplevels()
  
  # Sort by group_var to ensure consistent ordering
  df <- df %>% arrange(.data[[group_var]])
  
  # Extract trait values (Ly)
  Ly <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(vals = list(.data[[trait_col]]), .groups = "drop") %>%
    dplyr::pull(vals)
  
  # Extract time values (Lt)
  Lt <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(time = list(.data[[time_var]]), .groups = "drop") %>%
    dplyr::pull(time)
  
  # Run FPCA
  result <- fdapace::FPCA(Ly = Ly, Lt = Lt)
  
  # Store the group_var values used (in order) as an attribute
  group_ids <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data[[group_var]])
  
  attr(result, "group_ids") <- group_ids
  attr(result, "group_var") <- group_var
  
  return(result)
  
}

  

# Field data (blade)
fe_b_l <- run_fpca(my_dat1, group_var = "un_num", time_var = "Leaf", trait_col = "Length")   
fe_b_t <- run_fpca(my_dat1, group_var = "un_num", time_var = "rel_t", trait_col = "Length")
fe_b_e <- run_fpca(my_dat1, group_var = "un_num", time_var = "rel_e", trait_col = "Length")
fe_b_n <- run_fpca(my_dat1, group_var = "un_num", time_var = "rel_n", trait_col = "Length")

# ==============================================================================
# Shared helpers for sections 3a – 3c
# ==============================================================================

# ---- Extract FPC scores and join covariates from the original data ----
# origin_val : if non-NULL, the original_df is pre-filtered to that origin level
#              before computing the per-group t_i summary.
extract_fpca_scores <- function(fpca_result, original_df, group_var,
                                ti_var     = "t_i",
                                origin_val = NULL) {
  
  group_ids <- attr(fpca_result, "group_ids")
  scores    <- as.data.frame(fpca_result$xiEst)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores[[group_var]] <- group_ids
  
  ref_df <- if (!is.null(origin_val)) {
    dplyr::filter(original_df, origin == origin_val)
  } else {
    original_df
  }
  
  ti_df <- ref_df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      t_i    = mean(.data[[ti_var]], na.rm = TRUE),
      origin = dplyr::first(origin),
      .groups = "drop"
    )
  
  dplyr::left_join(scores, ti_df, by = group_var)
}

# ---- Regress each FPC score on t_i and return a tidy summary table ----
# Includes lambda_weight (λ_k / Σλ) so per-PC R² values can be meaningfully
# aggregated into a single composite R² without double-counting variance.
#
# origin_val : if non-NULL, scores are filtered to that origin before fitting.
#              This supports both the pooled-FPCA-split-by-origin approach (3b)
#              and the separate-FPCA-per-origin approach (3c).
regress_fpca_scores <- function(fpca_result, original_df, group_var,
                                ti_var     = "t_i",
                                label      = "",
                                origin_val = NULL) {
  
  scores_df <- extract_fpca_scores(fpca_result, original_df, group_var,
                                   ti_var, origin_val)
  
  if (!is.null(origin_val)) {
    scores_df <- dplyr::filter(scores_df, origin == origin_val)
  }
  
  pc_cols <- grep("^PC", names(scores_df), value = TRUE)
  lambdas <- fpca_result$lambda
  weights <- lambdas / sum(lambdas)
  
  purrr::map_dfr(seq_along(pc_cols), function(i) {
    pc      <- pc_cols[i]
    formula <- as.formula(paste(pc, "~ t_i"))
    fit     <- lm(formula, data = scores_df)
    fit_tidy <- broom::tidy(fit)
    
    broom::glance(fit) %>%
      dplyr::select(r.squared, adj.r.squared, statistic, p.value, df, df.residual) %>%
      dplyr::mutate(
        PC            = pc,
        analysis      = label,
        origin        = if (!is.null(origin_val)) origin_val else NA_character_,
        lambda_weight = weights[i],
        slope         = coef(fit)["t_i"],
        slope_p       = fit_tidy %>% dplyr::filter(term == "t_i") %>% dplyr::pull(p.value),
        .before       = 1
      )
  })
}

# ==============================================================================
# 3e. PC1 vs PC2 — linear predicted trajectory of T.I., by origin
#     PC1 and PC2 are each regressed on t_i separately within origin.
#     The joint predicted path (PC1_hat(t), PC2_hat(t)) is traced over a fine
#     T.I. grid, giving a continuous trajectory through the score space.
#     Confidence bands come from the prediction intervals of each model.
# ==============================================================================

# ---- Extract scores with t_i, reduce to genotype means per origin ----
# my_dat1 is long (one row per leaf), so pull the genotype lookup as distinct
# plant-level rows before joining, then average to genotype x origin.
geno_lookup <- my_dat1 %>%
  dplyr::distinct(un_num, genotype)

# ---- Helper: build + save one trajectory plot ----
plot_fpca_trajectory <- function(
    fe_model,          # one of fe_b_l, fe_b_t, fe_b_e, fe_b_n
    my_dat,            # data passed to extract_fpca_scores
    geno_lookup,       # lookup table with un_num -> variety_id_n, origin
    origin_colours,    # named colour vector
    fig_dir,
    timestamp,
    file_tag,          # short string for filename, e.g. "Fig5A"
    plot_title = "A",  # panel label
    # Optional per-origin, per-endpoint nudge overrides.
    # Pass a tibble with columns: origin, endpoint, nudge_x_mult, nudge_y_mult
    # These REPLACE (not add to) the defaults for matching rows.
    nudge_overrides = NULL,
    xlim = NULL,       # e.g. global_xlim
    ylim = NULL        # e.g. global_ylim
) {
  
  # ---- Build geno_df ----
  geno_df <- extract_fpca_scores(fe_model, my_dat, group_var = "un_num") %>%
    dplyr::filter(!is.na(t_i)) %>%
    dplyr::left_join(geno_lookup, by = "un_num") %>%
    dplyr::group_by(genotype, origin) %>%
    dplyr::summarise(
      PC1 = mean(PC1, na.rm = TRUE),
      PC2 = mean(PC2, na.rm = TRUE),
      t_i = mean(t_i, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ---- Fit linear trajectories per origin ----
  geno_origins <- unique(geno_df$origin)
  
  traj_df <- purrr::map_dfr(geno_origins, function(org) {
    df_org  <- dplyr::filter(geno_df, origin == org)
    fit_pc1 <- lm(PC1 ~ t_i, data = df_org)
    fit_pc2 <- lm(PC2 ~ t_i, data = df_org)
    
    ti_grid <- seq(min(df_org$t_i), max(df_org$t_i), length.out = 200)
    nd      <- data.frame(t_i = ti_grid)
    
    pred_pc1 <- predict(fit_pc1, newdata = nd, interval = "confidence") %>%
      as.data.frame() %>%
      dplyr::rename(PC1 = fit, PC1_lwr = lwr, PC1_upr = upr)
    
    pred_pc2 <- predict(fit_pc2, newdata = nd, interval = "confidence") %>%
      as.data.frame() %>%
      dplyr::rename(PC2 = fit, PC2_lwr = lwr, PC2_upr = upr)
    
    dplyr::bind_cols(nd, pred_pc1, pred_pc2) %>%
      dplyr::mutate(origin = org)
  })
  
  # ---- Endpoint nudge defaults (multipliers of axis range) ----
  pc1_range <- diff(range(geno_df$PC1))
  pc2_range <- diff(range(geno_df$PC2))
  
  # Default nudge table — edit multipliers here as a first pass;
  # use nudge_overrides for per-input fine-tuning.
  nudge_defaults <- tidyr::expand_grid(
    origin   = geno_origins,
    endpoint = c("early T.I.", "late T.I.")
  ) %>%
    dplyr::mutate(
      nudge_x_mult = dplyr::case_when(
        endpoint == "early T.I."  ~ -0.01,
        endpoint == "late T.I." ~  0.03,
        TRUE                    ~  0.0
      ),
      nudge_y_mult = dplyr::case_when(
        endpoint == "early T.I."  ~ -0.10,
        endpoint == "late T.I." ~  0.15,
        TRUE                    ~  0.0
      )
    )
  
  # Apply overrides: replace matching rows
  if (!is.null(nudge_overrides)) {
    nudge_defaults <- nudge_defaults %>%
      dplyr::rows_update(nudge_overrides, by = c("origin", "endpoint"))
  }
  
  endpoints <- traj_df %>%
    dplyr::group_by(origin) %>%
    dplyr::slice(c(1, dplyr::n())) %>%
    dplyr::mutate(endpoint = c("early T.I.", "late T.I.")) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(nudge_defaults, by = c("origin", "endpoint")) %>%
    dplyr::mutate(
      nudge_x = nudge_x_mult * pc1_range,
      nudge_y = nudge_y_mult * pc2_range
    )
  
  # ---- Arrow direction vectors ----
  arrow_df <- traj_df %>%
    dplyr::group_by(origin) %>%
    dplyr::slice_tail(n = 6) %>%
    dplyr::summarise(
      x    = first(PC1), y    = first(PC2),
      xend = last(PC1),  yend = last(PC2),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dx   = xend - x,   dy   = yend - y,
      xend = xend + 0.18 * dx,
      yend = yend + 0.18 * dy
    )
  
  # ---- Axis labels ----
  pve <- round(100 * fe_model$lambda / sum(fe_model$lambda), 0)
  panel_ratio <- pve[2] / pve[1]   # height / width ∝ FPC2 pve / FPC1 pve
  
  # ---- Plot ----
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(
      data  = geno_df,
      aes(x = PC1, y = PC2, colour = origin),
      size  = 1.8, alpha = 0.3
    ) +
    geom_ribbon(
      data  = traj_df,
      aes(x = PC1, ymin = PC2_lwr, ymax = PC2_upr, fill = origin),
      alpha = 0.2
    ) +
    geom_path(
      data = traj_df,
      aes(x = PC1, y = PC2, colour = origin),
      linewidth = 1.2, lineend = "round"
    ) +
    geom_segment(
      data = arrow_df,
      aes(x = x, y = y, xend = xend, yend = yend, colour = origin),
      linewidth = 1.2, lineend = "butt", linejoin = "mitre",
      arrow = arrow(type = "closed", length = unit(5, "mm"), angle = 16)
    ) +
    geom_text(
      data  = endpoints,
      aes(x = PC1 + nudge_x, y = PC2 + nudge_y, label = endpoint, colour = origin),
      size  = 5.0, show.legend = FALSE
    ) +
    scale_colour_manual(values = origin_colours, name = "Origin") +
    scale_fill_manual(values   = origin_colours, name = "Origin") +
    scale_y_continuous(expand  = expansion(mult = 0.1)) +
    labs(
      title = plot_title,
      x = paste0("FPC1 (", pve[1], "%)"),
      y = paste0("FPC2 (", pve[2], "%)")
    ) +
    guides(
      colour = guide_legend(override.aes = list(shape = c(16, 16), linetype = 0, size = 3)),
      fill   = "none"
    ) +
    # coord_equal(xlim = xlim, ylim = ylim, clip = "off") +
    # coord_cartesian(clip = "off", ratio = 1) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(
      aspect.ratio = panel_ratio,   # <-- toggle this for axes proportional to FPC PVE
      legend.position      = c(1.02, 1),
      legend.justification = c(0, 1),
      legend.title         = element_text(size = 11, face = "bold"),
      legend.text          = element_text(size = 12),
      plot.margin          = margin(10, 100, 10, 10),
      plot.title           = element_text(hjust = 0, size = 15, face = "bold"),
      axis.text.x          = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y          = element_text(size = 12),
      axis.title           = element_text(size = 14),
      axis.title.x         = element_text(size = 14, margin = margin(t = 10)),
      axis.title.y         = element_text(size = 14, margin = margin(r = 10)),
      panel.grid.minor     = element_line(linewidth = 0.3, colour = "grey85")
    )
  
  # ---- Save ----
  ggsave(
    filename = paste0(fig_dir, file_tag, "_", timestamp, ".pdf"),
    plot     = p,
    width    = 7.5, height = 5, dpi = 300
  )
  
  invisible(p)  # return silently so caller can inspect if needed
}

# ---- Run for all four inputs ----

# Define the four models with their metadata
model_specs <- list(
  list(model = fe_b_l, tag = "Fig5A", title = "A"),
  list(model = fe_b_t, tag = "Fig5B", title = "B"),
  list(model = fe_b_e, tag = "Fig5C", title = "C"),
  list(model = fe_b_n, tag = "Fig5D", title = "")
)

# Per-model nudge overrides — add rows here as you inspect each plot.
# NULL means "use defaults". Rows must match exactly on origin + endpoint.
nudge_list <- list(
  NULL,  # fe_b_l: defaults fine
  NULL,  # fe_b_t: defaults fine  
  # fe_b_e: example — push france late T.I. label further right
  tibble::tibble(
    origin = "france", endpoint = "late T.I.",
    nudge_x_mult = 0.06, nudge_y_mult = 0.10
  ),
  NULL   # fe_b_n: defaults fine
)

# ---- Colours: define once, used everywhere ----
# Adjust hex codes to taste
origin_colours <- c(
  "france" = "#0072B2",   # French blue
  "mexico" = "#D55E00"    # Mexican green
)


# ---- Run plots (individual, free-scale, saved) ----
plots <- purrr::map2(model_specs, nudge_list, function(spec, nudges) {
  plot_fpca_trajectory(
    fe_model        = spec$model,
    my_dat          = my_dat1,
    geno_lookup     = geno_lookup,
    origin_colours  = origin_colours,
    fig_dir         = fig_dir,
    timestamp       = timestamp,
    file_tag        = spec$tag,
    plot_title      = spec$title,
    nudge_overrides = nudges
    # no xlim/ylim — each panel uses its own natural scale
  )
})




