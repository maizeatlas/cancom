rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# T.I. functional variance explained: per-PC and composite R² across
# registration systems and origins (field experiment)
#
#   Pipeline:
#     1. FPCA is run on field blade-length profiles under four registration
#        systems (absolute rank, T.I.-relative, ear-relative, leaf-number
#        normalised); FPC scores are extracted per plant and joined to T.I.
#     2. Each FPC score is regressed on T.I. within origin (France, Mexico)
#        and pooled; per-PC adj. R² and slope significance are recorded
#        alongside eigenvalue weights (λk / Σλ).
#     3. A composite functional R² (Σk λk/Σλ · R²k) summarises the total
#        functional variance linearly attributable to T.I. timing, enabling
#        cross-system comparison without double-counting variance components.
#     4. Results are focused on the best-performing registration system
#        (Leaf No. Norm.) and visualised as lollipop plots of adj. R² per
#        FPC, faceted by subset (Pooled, France, Mexico), with composite R²
#        shown as a reference line per facet.
#
#   Outputs:
#     A. Fig panel A — lollipop chart of per-FPC adj. R² for Leaf No. Norm.,
#        faceted by subset, with significance stars and composite R² overlay;
#        saved as Fig_compR2_lollipop.pdf.
# ==============================================================================

library(tidyverse) #data manipulation
library(broom) #model output
library(fdapace) #functional data analysis
library(patchwork) #combining plots


# ---- Paths ----
dat_dir   <- paste0(getwd(), "/data/")
fig_dir   <- paste0(getwd(), "/figures/")
timestamp <- format(Sys.time(), "%d%m%Y")

# ==============================================================================
# 1. Read data
# ==============================================================================

#canopy dimensions from the field
name_fixes <- read_xlsx(paste0(dat_dir, "pg_unmatched_gtypes_in_ptypes.xlsx")) %>%
  dplyr::mutate(sample_ptype = tolower(sample_ptype))

my_dat1 <- read.csv(paste0(dat_dir, "05_dat_field_canopy.csv")) %>%
  #dplyr::filter(un_num != "98_4_1" & un_num != "122_1_2")
  dplyr::rename(genotype = variety_id_n) %>%
  dplyr::mutate(genotype = tolower(genotype)) %>%
  dplyr::left_join(name_fixes, by = c("genotype" = "sample_ptype")) %>%
  dplyr::mutate(genotype = dplyr::coalesce(sample_gtype, genotype)) %>%
  dplyr::select(-sample_gtype) %>%
  dplyr::filter(genotype != "NA")

#quick check
sfd<-my_dat1%>%
  dplyr::filter(!is.na(Length))%>%
  droplevels()%>%
  # dplyr::group_by(un_num)%>%
  # dplyr::select(un_num)%>%
  # distinct()%>%
  #  group_by(un_num)%>%
  #   tally()
  dplyr::group_by(origin)%>%
  dplyr::select(variety_id_n)%>%
  distinct()%>%
  count()



# ==============================================================================
# 2. FPCA
# ==============================================================================

run_fpca <- function(df, group_var, time_var, trait_col) {
  df <- df %>%
    dplyr::filter(!is.na(.data[[trait_col]]), !is.na(.data[[time_var]])) %>%
    droplevels() %>%
    dplyr::arrange(.data[[group_var]])
  
  Ly <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(vals = list(.data[[trait_col]]), .groups = "drop") %>%
    dplyr::pull(vals)
  
  Lt <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(time = list(.data[[time_var]]), .groups = "drop") %>%
    dplyr::pull(time)
  
  result    <- fdapace::FPCA(Ly = Ly, Lt = Lt)
  group_ids <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data[[group_var]])
  
  attr(result, "group_ids") <- group_ids
  attr(result, "group_var") <- group_var
  result
}

fe_b_l <- run_fpca(my_dat1, "un_num", "Leaf",  "Length")
fe_b_t <- run_fpca(my_dat1, "un_num", "rel_t", "Length")
fe_b_e <- run_fpca(my_dat1, "un_num", "rel_e", "Length")
fe_b_n <- run_fpca(my_dat1, "un_num", "rel_n", "Length")

# ==============================================================================
# 3. Helpers
# ==============================================================================

extract_fpca_scores <- function(fpca_result, original_df, group_var,
                                ti_var = "t_i", origin_val = NULL) {
  group_ids <- attr(fpca_result, "group_ids")
  scores    <- as.data.frame(fpca_result$xiEst)
  colnames(scores) <- paste0("PC", seq_len(ncol(scores)))
  scores[[group_var]] <- group_ids
  
  ref_df <- if (!is.null(origin_val)) dplyr::filter(original_df, origin == origin_val) else original_df
  
  ti_df <- ref_df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      t_i    = mean(.data[[ti_var]], na.rm = TRUE),
      origin = dplyr::first(origin),
      .groups = "drop"
    )
  
  dplyr::left_join(scores, ti_df, by = group_var)
}

regress_fpca_scores <- function(fpca_result, original_df, group_var,
                                ti_var = "t_i", label = "", origin_val = NULL) {
  scores_df <- extract_fpca_scores(fpca_result, original_df, group_var,
                                   ti_var, origin_val)
  if (!is.null(origin_val)) scores_df <- dplyr::filter(scores_df, origin == origin_val)
  
  pc_cols <- grep("^PC", names(scores_df), value = TRUE)
  weights <- fpca_result$lambda / sum(fpca_result$lambda)
  
  purrr::map_dfr(seq_along(pc_cols), function(i) {
    pc      <- pc_cols[i]
    fit     <- lm(as.formula(paste(pc, "~ t_i")), data = scores_df)
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

summarise_functional_r2 <- function(results_df) {
  results_df %>%
    dplyr::group_by(analysis, origin) %>%
    dplyr::summarise(functional_r2 = sum(lambda_weight * r.squared), .groups = "drop")
}

sig_label <- function(p) {
  dplyr::case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

# ==============================================================================
# 4. Regressions: pooled and split by origin
# ==============================================================================

analysis_order <- c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm.")
origins        <- unique(my_dat1$origin)

fe_pooled <- dplyr::bind_rows(
  regress_fpca_scores(fe_b_l, my_dat1, "un_num", label = "Abs. Rank"),
  regress_fpca_scores(fe_b_t, my_dat1, "un_num", label = "Rel. T.I."),
  regress_fpca_scores(fe_b_e, my_dat1, "un_num", label = "Rel. Ear Rank"),
  regress_fpca_scores(fe_b_n, my_dat1, "un_num", label = "Leaf No. Norm.")
) %>%
  dplyr::mutate(
    analysis = factor(analysis, levels = analysis_order),
    PC       = forcats::fct_relabel(PC, ~ sub("^PC", "FPC", .x))
  )

fe_pooled_compR2 <- summarise_functional_r2(fe_pooled)

fe_pooled_splitby_origin <- dplyr::bind_rows(
  purrr::map_dfr(origins, ~ regress_fpca_scores(fe_b_l, my_dat1, "un_num", label = "Abs. Rank",      origin_val = .x)),
  purrr::map_dfr(origins, ~ regress_fpca_scores(fe_b_t, my_dat1, "un_num", label = "Rel. T.I.",      origin_val = .x)),
  purrr::map_dfr(origins, ~ regress_fpca_scores(fe_b_e, my_dat1, "un_num", label = "Rel. Ear Rank",  origin_val = .x)),
  purrr::map_dfr(origins, ~ regress_fpca_scores(fe_b_n, my_dat1, "un_num", label = "Leaf No. Norm.", origin_val = .x))
) %>%
  dplyr::mutate(
    analysis = factor(analysis, levels = analysis_order),
    PC       = forcats::fct_relabel(PC, ~ sub("^PC", "FPC", .x))
  )

fe_pooled_splitby_origin_comR2 <- summarise_functional_r2(fe_pooled_splitby_origin)

# ==============================================================================
# 5. Lollipop plot — Leaf No. Norm. focused panel
# ==============================================================================

target_method <- "Leaf No. Norm."

ln_pc <- dplyr::bind_rows(
  fe_pooled %>%
    dplyr::filter(analysis == target_method) %>%
    dplyr::mutate(subset = "Pooled"),
  fe_pooled_splitby_origin %>%
    dplyr::filter(analysis == target_method) %>%
    dplyr::mutate(subset = tools::toTitleCase(as.character(origin)))
) %>%
  dplyr::transmute(subset, PC, r2 = adj.r.squared, sig = sig_label(slope_p)) %>%
  dplyr::mutate(
    subset = factor(subset, levels = c("Pooled", "France", "Mexico")),
    PC     = factor(PC, levels = rev(unique(PC)))
  )

ln_comp <- dplyr::bind_rows(
  fe_pooled_compR2 %>%
    dplyr::filter(analysis == target_method) %>%
    dplyr::mutate(subset = "Pooled"),
  fe_pooled_splitby_origin_comR2 %>%
    dplyr::filter(analysis == target_method) %>%
    dplyr::mutate(subset = tools::toTitleCase(as.character(origin)))
) %>%
  dplyr::transmute(
    subset = factor(subset, levels = c("Pooled", "France", "Mexico")),
    functional_r2
  )

subset_cols <- c(Pooled = "grey45", France = "#1f77b4", Mexico = "#d62728")

ln_focus_plot <- ggplot(ln_pc, aes(x = r2, y = PC, color = subset)) +
  geom_vline(data = ln_comp, aes(xintercept = functional_r2),
             linetype = "dashed", color = "grey55", linewidth = 0.4) +
  geom_text(data = ln_comp,
            aes(x = functional_r2,
                y = levels(ln_pc$PC)[length(levels(ln_pc$PC))],
                label = sprintf("italic(R)[c]^2 == %.2f", functional_r2)),
            parse = TRUE, inherit.aes = FALSE,
            size = 3.25, hjust = -0.05, vjust = -0.6, color = "grey30") +
  geom_segment(aes(x = 0, xend = r2, yend = PC), linewidth = 0.5) +
  geom_point(size = 2) +
  geom_text(aes(label = sig, vjust = ifelse(sig == "ns", 0.25, 0.75)),
            nudge_x = 0.1, size = 3, color = "grey15") +
  facet_wrap(~ subset, nrow = 1) +
  scale_color_manual(values = subset_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 1.15), breaks = c(0, 0.5, 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = c(0.5, 1.8))) +
  labs(title = "A", x = expression(Adj-R^2), y = NULL) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey92", color = NA),
    strip.text         = element_text(face = "bold", size = 10),
    plot.title         = element_text(face = "bold", size = 12),
    axis.text.x        = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y        = element_text(size = 10),
    axis.title         = element_text(size = 12),
    axis.title.x       = element_text(size = 12, margin = margin(t = 6)),
    axis.title.y       = element_text(size = 12, margin = margin(r = 6)),
    plot.margin        = margin(3, 3, 3, 3)
  )

ln_focus_plot

ggsave(paste0(fig_dir, "Fig_compR2_lollipop_", timestamp, ".pdf"),
       ln_focus_plot,
       width = 5.5, height = 2, dpi = 300)
