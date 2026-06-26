rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# Duration-independent leaf scaling: T.I.-corrected longest-leaf residuals
#
#   Pipeline:
#     1. FPCA reconstructs blade-length profiles per plant (rel_n time system).
#     2. Longest leaf is extracted per plant and averaged to genotype x env.
#     3. The residual of (longest leaf ~ T.I.), fit WITHIN each environment,
#        gives the scaling variation NOT attributable to developmental duration
#        (i.e. to T.I. timing).
#
#   Outputs:
#     A. GWAS input file of residual scaling values (per env + combined mean).
#     B. France-vs-Mexico residual scatter classifying genotypes as
#        hyper- / hypo-scalers (consistent across both envs) or env-specific.
# ==============================================================================

library(tidyverse)   # dplyr, tidyr, ggplot2, purrr, tibble, stringr
library(fdapace)     # FPCA
library(readxl)      # read_xlsx (genotype name fixes)

# ---- Paths ----
dat_dir   <- paste0(getwd(), "/data/")
fig_dir   <- paste0(getwd(), "/figures/")
timestamp <- format(Sys.time(), "%d%m%Y")

# ==============================================================================
# 1. Read data & harmonise genotype names
# ==============================================================================
name_fixes <- read_xlsx(paste0(dat_dir, "pg_unmatched_gtypes_in_ptypes.xlsx")) %>%
  dplyr::mutate(sample_ptype = tolower(sample_ptype))

my_dat1 <- read.csv(paste0(dat_dir, "05_dat_field_canopy.csv")) %>%
  dplyr::rename(genotype = variety_id_n) %>%
  dplyr::mutate(genotype = tolower(genotype)) %>%
  dplyr::left_join(name_fixes, by = c("genotype" = "sample_ptype")) %>%
  dplyr::mutate(genotype = dplyr::coalesce(sample_gtype, genotype)) %>%
  dplyr::select(-sample_gtype) %>%
  dplyr::filter(genotype != "NA")

# ==============================================================================
# 2. FPCA workflow (reconstruct blade-length profiles)
# ==============================================================================
run_fpca_workflow <- function(df, group_var, time_var, trait_col,
                              fpca_optns = list(FVEthreshold = 0.999), ...) {
  
  # Reshape to one (Ly, Lt) list per subject
  input_tbl <- df %>%
    dplyr::filter(!is.na(.data[[trait_col]]), !is.na(.data[[time_var]])) %>%
    dplyr::arrange(.data[[group_var]], .data[[time_var]]) %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      Ly = list(.data[[trait_col]]),
      Lt = list(.data[[time_var]]),
      .groups = "drop"
    ) %>%
    dplyr::rename(subject = dplyr::all_of(group_var))
  
  Ly <- input_tbl$Ly; Lt <- input_tbl$Lt; subject_ids <- input_tbl$subject
  
  fit <- fdapace::FPCA(Ly = Ly, Lt = Lt, optns = fpca_optns, ...)
  
  fitted_mat <- {
    fo <- fitted(fit)
    if (is.list(fo) && !is.null(fo$fitted)) fo$fitted else fo
  }
  
  # Reconstructed curves on the common work grid
  recon_df <- purrr::map_dfr(seq_along(Ly), function(i) {
    tibble::tibble(subject = subject_ids[i], time = fit$workGrid, recon = fitted_mat[i, ])
  })
  
  list(fit = fit, recon_df = recon_df)
}

fpca_out <- run_fpca_workflow(
  df = my_dat1, group_var = "un_num", time_var = "rel_n", trait_col = "Length"
)

# ==============================================================================
# 3. Longest leaf per plant -> genotype x env means
#    Curves are mapped from the leaf-normalised grid back to leaf-rank units
#    (time * LN_total), then the ascending early ranks (< 4.5) are dropped
#    before taking the maximum blade length.
# ==============================================================================
meta <- my_dat1 %>%
  dplyr::select(un_num, genotype, origin, LN_total, t_i) %>%
  dplyr::distinct() %>%
  dplyr::rename(subject = un_num)

longest_leaf <- fpca_out$recon_df %>%
  dplyr::left_join(meta, by = "subject") %>%
  dplyr::mutate(time = LN_total * time) %>%
  dplyr::filter(time >= 4.5) %>%
  dplyr::group_by(subject, genotype, origin, t_i) %>%
  dplyr::summarise(longest_leaf = max(recon, na.rm = TRUE), .groups = "drop")

# Genotype x env means (env capitalised to France / Mexico)
geno_means <- longest_leaf %>%
  dplyr::rename(env = origin) %>%
  dplyr::group_by(genotype, env) %>%
  dplyr::summarise(
    meanti     = mean(t_i,         na.rm = TRUE),
    meanlength = mean(longest_leaf, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(env = str_to_sentence(env))

# Fit per-env models and build equation labels (one per facet)
eqn_labels <- geno_means %>%
  dplyr::group_by(env) %>%
  dplyr::group_modify(~ {
    m  <- lm(meanlength ~ meanti, data = .x)
    cf <- coef(m)
    r2 <- summary(m)$r.squared
    p  <- summary(m)$coefficients[2, 4]
    stars <- dplyr::case_when(
      p < 0.001 ~ "'***'", p < 0.01 ~ "'**'", p < 0.05 ~ "'*'", TRUE ~ "'ns'"
    )
    tibble::tibble(
      label = sprintf("italic(y) == %.2f + %.2f * italic(x)", cf[1], cf[2]),
      r2lab = sprintf("italic(R)^2 == %.2f~~%s", r2, stars)
    )
  }) %>%
  dplyr::ungroup()

env_cols <- c(France = "#0072B2", Mexico = "#D55E00")

p_ti_length <- ggplot(geno_means, aes(meanti, meanlength, color = env)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
  # equation (top line) and R2 + significance (second line)
  geom_text(data = eqn_labels, aes(x = -Inf, y = Inf, label = label),
            parse = TRUE, hjust = -0.05, vjust = 1.5, size = 3.2,
            color = "grey15", inherit.aes = FALSE) +
  geom_text(data = eqn_labels, aes(x = -Inf, y = Inf, label = r2lab),
            parse = TRUE, hjust = -0.05, vjust = 3.2, size = 3.2,
            color = "grey15", inherit.aes = FALSE) +
  facet_wrap(~ env) +
  scale_color_manual(values = env_cols, guide = "none") +
  labs(
    x = "T.I. timing",
    y = "Longest leaf length (cm)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    axis.title.x     = element_text(margin = margin(t = 6)),
    axis.title.y     = element_text(margin = margin(r = 6))
  )
p_ti_length

# ggsave(paste0(fig_dir, "Fig_ti_length_regression_", timestamp, ".pdf"),
#        p_ti_length, width = 7, height = 3.5, dpi = 300)

# Compute residuals per environment
resid_dat <- geno_means %>%
  dplyr::group_by(env) %>%
  dplyr::mutate(
    fitted = predict(lm(meanlength ~ meanti)),
    resid  = residuals(lm(meanlength ~ meanti))
  ) %>%
  dplyr::ungroup()

# Optional: test for heteroscedasticity (Breusch-Pagan) per env, for annotation
bp_labels <- resid_dat %>%
  dplyr::group_by(env) %>%
  dplyr::group_modify(~ {
    m  <- lm(meanlength ~ meanti, data = .x)
    bp <- lmtest::bptest(m)               # requires lmtest
    tibble::tibble(label = sprintf("BP~italic(p) == %.2f", bp$p.value))
  }) %>%
  dplyr::ungroup()

env_cols <- c(France = "#0072B2", Mexico = "#D55E00")

p_resid <- ggplot(resid_dat, aes(meanti, resid, color = env)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.7) +
  # smoother to reveal any trend in spread/level (heteroscedasticity shows as curvature/fanning)
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.7, alpha = 0.15, color = "grey30") +
  geom_text(data = bp_labels, aes(x = -Inf, y = Inf, label = label),
            parse = TRUE, hjust = -0.1, vjust = 1.5, size = 3,
            color = "grey15", inherit.aes = FALSE) +
  facet_wrap(~ env) +
  scale_color_manual(values = env_cols, guide = "none") +
  labs(
    x = "T.I. timing",
    y = "Residual longest leaf length (cm) | T.I."
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 10),
    axis.title.x     = element_text(margin = margin(t = 6)),
    axis.title.y     = element_text(margin = margin(r = 6))
  )
p_resid

# ggsave(paste0(fig_dir, "Fig_residual_ti_", timestamp, ".pdf"),
#        p_resid, width = 7, height = 3.5, dpi = 300)

# ==============================================================================
# 4. T.I. correction: residual of longest leaf ~ T.I., fit WITHIN each env
#    + residual (longer than predicted by T.I.) -> over-scaling  ("hyper")
#    - residual (shorter than predicted)        -> under-scaling ("hypo")
# ==============================================================================

# 4a. Genotype-mean level (one residual per genotype × env)
geno_resid <- geno_means %>%
  dplyr::group_by(env) %>%
  dplyr::mutate(resid_length_ti = residuals(lm(meanlength ~ meanti))) %>%
  dplyr::ungroup()

# 4b. Per-plant level (one residual per plant; plant-level within-env fit)
plant_resid <- longest_leaf %>%
  dplyr::rename(env = origin) %>%
  dplyr::mutate(env = str_to_sentence(env)) %>%
  dplyr::group_by(env) %>%
  dplyr::mutate(resid_length_ti = residuals(lm(longest_leaf ~ t_i))) %>%
  dplyr::ungroup()

# ==============================================================================
# 5. OUTPUT A — GWAS residual files (genotype means + per-plant)
# ==============================================================================

# 5a. Genotype-mean residuals (per-env rows + combined mean across envs)
resid_table <- dplyr::bind_rows(
  geno_resid %>%
    dplyr::transmute(genotype, origin = tolower(env), value = resid_length_ti),
  geno_resid %>%
    dplyr::group_by(genotype) %>%
    dplyr::summarise(value = mean(resid_length_ti, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(origin = "combined")
) %>%
  dplyr::transmute(
    genotype, origin,
    trait_name  = "resid_longlength_ti",
    value,
    fpca_source = NA_character_,
    trait       = "mean"
  ) %>%
  dplyr::arrange(genotype, origin)

write.table(resid_table, paste0(dat_dir, "06_dat_gwas_scaling_residvar.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 5b. Per-plant residuals (subject = plot_plant_env, e.g. "100_1_2")
plant_resid_table <- plant_resid %>%
  dplyr::transmute(
    subject, genotype,
    origin = tolower(env),
    value  = resid_length_ti,
    trait_name  = "resid_longlength_ti_plant",
    fpca_source = NA_character_,
    trait       = "individual"
  ) %>%
  tidyr::separate(
    subject,
    into   = c("plot", "plant", "env_code"),
    sep    = "_",
    remove = FALSE          # keep the original subject column
  ) %>%
  dplyr::arrange(genotype, origin, plot, plant)

write.table(plant_resid_table,
            paste0(dat_dir, "06_dat_gwas_scaling_residvar_perplant.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# 6. Variation summary: spread before vs after T.I. correction, per env
#    var_reduction_pct equals the per-env R^2 of lm(meanlength ~ meanti).
# ==============================================================================
var_compare <- geno_resid %>%
  dplyr::group_by(origin = env) %>%
  dplyr::summarise(
    var_before   = var(meanlength,             na.rm = TRUE),
    var_after    = var(resid_length_ti,        na.rm = TRUE),
    sd_before    = sd(meanlength,              na.rm = TRUE),
    sd_after     = sd(resid_length_ti,         na.rm = TRUE),
    range_before = diff(range(meanlength,      na.rm = TRUE)),
    range_after  = diff(range(resid_length_ti, na.rm = TRUE)),
    n            = sum(!is.na(meanlength)),
    .groups      = "drop"
  ) %>%
  dplyr::mutate(
    var_reduction_pct   = 100 * (var_before   - var_after)   / var_before,
    range_reduction_pct = 100 * (range_before - range_after) / range_before
  )
print(var_compare)

# ==============================================================================
# 7. Wide residuals + quadrant classification (single canonical table)
#    FLIP the two diagonal labels in case_when if your sign convention differs.
# ==============================================================================
resid_wide <- geno_resid %>%
  dplyr::select(genotype, env, resid_length_ti) %>%
  tidyr::pivot_wider(names_from = env, values_from = resid_length_ti) %>%
  dplyr::filter(!is.na(France), !is.na(Mexico)) %>%   # genotypes present in both
  dplyr::mutate(
    quadrant = dplyr::case_when(
      France > 0 & Mexico > 0 ~ "Hyper (both)",
      France < 0 & Mexico < 0 ~ "Hypo (both)",
      TRUE                    ~ "Env-specific"
    ),
    mean_resid = (France + Mexico) / 2
  ) %>%
  dplyr::arrange(dplyr::desc(mean_resid))

r_cross <- cor(resid_wide$France, resid_wide$Mexico, use = "complete.obs")

# ==============================================================================
# 8. Scaler genotype lists (for downstream / other use)
# ==============================================================================
# 8a. BETWEEN-env (consistent): same sign in France AND Mexico
scalers_between <- resid_wide %>%
  dplyr::filter(quadrant %in% c("Hyper (both)", "Hypo (both)")) %>%
  dplyr::transmute(
    genotype,
    scaler_class = dplyr::recode(quadrant, "Hyper (both)" = "hyper", "Hypo (both)" = "hypo"),
    France, Mexico, mean_resid
  ) %>%
  dplyr::arrange(dplyr::desc(mean_resid))

# 8b. WITHIN-env: per genotype, per env, sign-based class (long format)
scalers_within <- geno_resid %>%
  dplyr::filter(genotype %in% resid_wide$genotype) %>%
  dplyr::transmute(
    genotype, env, resid_length_ti,
    scaler_class = dplyr::if_else(resid_length_ti > 0, "hyper", "hypo")
  ) %>%
  dplyr::arrange(env, dplyr::desc(resid_length_ti))

# Console summary
cat(sprintf(
  "Cross-env residual correlation r = %.2f\nGenotypes in both envs: %d | Hyper(both): %d | Hypo(both): %d | Env-specific: %d\n",
  r_cross, nrow(resid_wide),
  sum(resid_wide$quadrant == "Hyper (both)"),
  sum(resid_wide$quadrant == "Hypo (both)"),
  sum(resid_wide$quadrant == "Env-specific")
))

scaler_counts <- resid_wide %>%
  dplyr::count(quadrant, name = "n") %>%
  dplyr::mutate(pct = round(100 * n / sum(n), 1))
scaler_counts

# Strongest few at each consistent end (e.g. for labelling another figure)
top_hyper <- dplyr::slice_max(dplyr::filter(scalers_between, scaler_class == "hyper"), mean_resid, n = 5)
top_hypo  <- dplyr::slice_min(dplyr::filter(scalers_between, scaler_class == "hypo"),  mean_resid, n = 5)

# ==============================================================================
# 9. OUTPUT B — France vs Mexico residual scatter (hyper / hypo scaler classes)
# ==============================================================================
hyper_col <- "#2166AC"   # over-scaling  (upper-right quadrant)
hypo_col  <- "#B2182B"   # under-scaling (lower-left quadrant)
gxe_col   <- "grey55"    # environment-specific (off-diagonal)
quad_cols <- c("Hyper (both)" = hyper_col, "Hypo (both)" = hypo_col, "Env-specific" = gxe_col)

base_sz   <- 10
lim       <- max(abs(c(resid_wide$France, resid_wide$Mexico)), na.rm = TRUE) * 1.10
lab_inset <- lim * 0.97   # corner-label placement, just inside the panel

scaler_scatter <- ggplot(resid_wide, aes(France, Mexico)) +
  # consistent-quadrant shading (diagonal corners)
  annotate("rect", xmin = 0,    xmax = Inf, ymin = 0,    ymax = Inf, fill = hyper_col, alpha = 0.07) +
  annotate("rect", xmin = -Inf, xmax = 0,   ymin = -Inf, ymax = 0,   fill = hypo_col,  alpha = 0.07) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = quadrant), size = 2, alpha = 0.85) +
  # corner annotations: 2 consistent (diagonal) + 2 env-specific (off-diagonal)
  annotate("text", x =  lab_inset, y =  lab_inset, hjust = 1, vjust = 1,
           label = "Hyper-scalers\n(both env)", color = hyper_col,
           fontface = "bold", size = 2.9, lineheight = 0.9) +
  annotate("text", x = -lab_inset, y = -lab_inset, hjust = 0, vjust = 0,
           label = "Hypo-scalers\n(both env)", color = hypo_col,
           fontface = "bold", size = 2.9, lineheight = 0.9) +
  annotate("text", x = -lab_inset, y =  lab_inset, hjust = 0, vjust = 1,
           label = "Hypo FR /\nHyper MX", color = "grey45",
           fontface = "italic", size = 2.9, lineheight = 0.9) +
  annotate("text", x =  lab_inset, y = -lab_inset, hjust = 1, vjust = 0,
           label = "Hyper FR /\nHypo MX", color = "grey45",
           fontface = "italic", size = 2.9, lineheight = 0.9) +
  scale_color_manual(values = quad_cols, name = NULL) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    title = "",
    x = "Longest blade length (cm) | T.I., France",
    y = "Longest blade length (cm) | T.I., Mexico"
  ) +
  theme_bw(base_size = base_sz) +
  theme(
    aspect.ratio       = 1,
    panel.grid.minor   = element_line(linewidth = 0.3, colour = "grey85"),
    plot.title         = element_text(hjust = 0, vjust = -5, size = 12, face = "bold"),
    legend.title       = element_text(size = base_sz, face = "bold"),
    legend.text        = element_text(size = base_sz),
    legend.position    = "top",
    legend.box.spacing = unit(5, "pt"),
    legend.margin      = margin(0, 0, 0, 0),
    axis.text          = element_text(size = base_sz),
    axis.title         = element_text(size = base_sz + 2),
    axis.title.x       = element_text(size = base_sz + 2, margin = margin(t = 10)),
    axis.title.y       = element_text(size = base_sz + 2, margin = margin(r = 10)),
    plot.margin        = margin(10, 10, 10, 10)
  )
scaler_scatter

ggsave(paste0(fig_dir, "Fig6_", timestamp, ".pdf"),
       scaler_scatter, width = 5, height = 5.4, dpi = 300)



