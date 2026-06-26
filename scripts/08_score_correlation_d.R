rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# FPCA feature extraction: reconstructed blade profiles to canopy descriptors
# under competing registration systems (field experiment)
#
#   Pipeline:
#     1. A generalised FPCA workflow (run_fpca_workflow) fits fdapace per plant,
#        returning FPC scores, per-observation residuals, and reconstructed
#        curves on the work grid; a wrapper (run_fpca_systems) applies this
#        across all four registration systems in a single call.
#     2. Registration system fit quality is benchmarked via compare_fpca_systems,
#        which computes total absolute residuals per plant and residual-by-time
#        profiles across systems, printing a mean ± SD summary table.
#     3. From the best-fitting system (rel_n), extract_fpca_features derives
#        four sets of outputs per plant: FPC scores, reconstructed curves with
#        metadata, canopy descriptors (longest leaf length and position, scaling
#        rate, centre-of-mass position along the phytomer axis), and ascending-
#        limb derivatives (dy/dx vs. leaf rank) up to the longest-leaf peak.
#     4. score_cor_table correlates each FPC score against the canopy
#        descriptors (longest leaf, longest-leaf position, centre of mass) via
#        Pearson tests, returning r, p, and significance stars per FPC × trait
#        combination.
#
#   Outputs:
#     A. features$scores       — FPC scores per plant under rel_n registration.
#     B. features$recon        — reconstructed blade curves with plant metadata.
#     C. features$longest_leaf — per-plant canopy descriptors (length, position,
#                                scaling rate, centre of mass).
#     D. features$derivs       — ascending-limb slope profiles (feeds downstream
#                                max-slope window analysis).
#     E. score_cor_table()     — FPC × canopy-descriptor correlation table.
# ==============================================================================

library(tidyverse)#data manipulation
library(data.table)#data manipulation
library(zoo)#data manipulation
library(ggh4x) #facetting
library(patchwork)#combining plots
library(cowplot) #combining plots
library(fdapace) #functional data analysis
library(ggpmisc) #plot stats


# ---- Paths ----
dat_dir= paste0(getwd(),"/data/")
dat_dir_r = paste0(getwd(),"/data/raw_data/")
fig_dir = paste0(getwd(),"/figures/")
timestamp <- format(Sys.time(), "%d%m%Y")

# ==============================================================================
# 1. Read and summarize data ----------------
# ==============================================================================


name_fixes <- read_xlsx(paste0(dat_dir,"pg_unmatched_gtypes_in_ptypes.xlsx")) %>%
  mutate(sample_ptype = tolower(sample_ptype))


my_dat1 <- read.csv(paste0(dat_dir,"05_dat_field_canopy.csv")) %>%
  dplyr::rename(genotype = variety_id_n) %>% 
  mutate(genotype = tolower(genotype)) %>% 
  left_join(name_fixes, by = c("genotype" = "sample_ptype")) %>%
  mutate(genotype = coalesce(sample_gtype, genotype)) %>%
  dplyr::select(-sample_gtype) %>%
  dplyr::filter(genotype != "NA") #%>%
# dplyr::filter(origin=="france")



# ==============================================================================
# 2. Reonstruct curves----------------
# ==============================================================================
run_fpca_workflow <- function(df,
                              group_var,
                              time_var,
                              trait_col,
                              ti_col      = NULL,
                              fpca_optns  = list(FVEthreshold = 0.999),
                              ...) {
  
  # ── 1. Prepare input ────────────────────────────────────────────────────────
  df2 <- df %>%
    dplyr::filter(
      !is.na(.data[[trait_col]]),
      !is.na(.data[[time_var]])
    )
  
  if (!is.null(ti_col)) {
    df2 <- df2 %>% dplyr::filter(!is.na(.data[[ti_col]]))
  }
  
  input_tbl <- df2 %>%
    dplyr::arrange(.data[[group_var]], .data[[time_var]]) %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      Ly = list(.data[[trait_col]]),
      Lt = list(.data[[time_var]]),
      .groups = "drop"
    ) %>%
    dplyr::rename(subject = dplyr::all_of(group_var))
  
  Ly          <- input_tbl$Ly
  Lt          <- input_tbl$Lt
  subject_ids <- input_tbl$subject
  
  # ── 2. Fit FPCA ─────────────────────────────────────────────────────────────
  fit <- fdapace::FPCA(Ly = Ly, Lt = Lt, optns = fpca_optns, ...)
  
  fitted_mat <- {
    fo <- fitted(fit)
    if (is.list(fo) && !is.null(fo$fitted)) fo$fitted else fo
  }
  
  # ── 3. Residuals at observed points ─────────────────────────────────────────
  errs_df <- purrr::map_dfr(seq_along(Ly), function(i) {
    pred <- approx(fit$workGrid, fitted_mat[i, ], xout = Lt[[i]], rule = 2)$y
    tibble::tibble(
      subject = subject_ids[i],
      time    = Lt[[i]],
      obs     = Ly[[i]],
      pred    = pred,
      resid   = Ly[[i]] - pred,
      n_obs   = length(Ly[[i]])
    )
  })
  
  # ── 4. Scores ────────────────────────────────────────────────────────────────
  scores_df <- as.data.frame(fit$xiEst) %>%
    setNames(paste0("fpc", seq_len(ncol(.)))) %>%
    tibble::add_column(subject = subject_ids, .before = 1)
  
  # ── 5. Reconstructed curves on workGrid ─────────────────────────────────────
  recon_df <- purrr::map_dfr(seq_along(Ly), function(i) {
    tibble::tibble(
      subject = subject_ids[i],
      time    = fit$workGrid,
      recon   = fitted_mat[i, ]
    )
  })
  
  list(
    fit         = fit,
    subject_ids = subject_ids,
    Ly          = Ly,
    Lt          = Lt,
    errs_df     = errs_df,
    scores_df   = scores_df,
    recon_df    = recon_df
  )
}

compare_fpca_systems <- function(df,
                                 group_var,
                                 time_vars,
                                 trait_col,
                                 ti_col = NULL,
                                 ...) {
  
  results <- lapply(names(time_vars), function(sys) {
    cat(sprintf("\nRunning FPCA for system: %s\n", sys))
    run_fpca_workflow(
      df        = df,
      group_var = group_var,
      time_var  = time_vars[[sys]],
      trait_col = trait_col,
      ti_col    = ti_col,
      ...
    )
  })
  names(results) <- names(time_vars)
  
  # ── 2. Total absolute residual per subject per system (panel A) ───────────────
  error_total <- lapply(names(results), function(sys) {
    results[[sys]]$errs_df %>%
      dplyr::mutate(system = sys) %>%
      dplyr::group_by(subject, system) %>%
      dplyr::summarise(
        tot_resid_s = sum(abs(resid), na.rm = TRUE),
        .groups = "drop"
      )
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(system = dplyr::recode(system,
                                         "abs_rnk" = "Abs. Rank",
                                         "rel_ti"  = "Rel. T.I.",
                                         "rel_e"   = "Rel. Ear Rank",
                                         "rel_n"   = "Leaf No. Norm."
    )) %>%
    dplyr::mutate(system = factor(system,
                                  levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm.")))
  
  # ── 3. Residual by time per system (panel B) ──────────────────────────────────
  resid_time <- lapply(names(results), function(sys) {
    results[[sys]]$errs_df %>%
      dplyr::mutate(system = sys)
    #dplyr::group_by(system, time) %>%
    #dplyr::summarise(
    #  mean_resid = mean(resid, na.rm = TRUE),
    #  n          = dplyr::n(),
    #  .groups    = "drop"
    #)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(system = dplyr::recode(system,
                                         "abs_rnk" = "Abs. Rank",
                                         "rel_ti"  = "Rel. T.I.",
                                         "rel_e"   = "Rel. Ear Rank",
                                         "rel_n"   = "Leaf No. Norm."
    )) %>%
    dplyr::mutate(system = factor(system,
                                  levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm.")))
  
  # ── 4. Print summary ──────────────────────────────────────────────────────────
  cat("\n── Total absolute residual by system ──\n")
  error_total %>%
    dplyr::group_by(system) %>%
    dplyr::summarise(
      mean = mean(tot_resid_s, na.rm = TRUE),
      sd   = sd(tot_resid_s,   na.rm = TRUE),
      .groups = "drop"
    ) %>%
    print()
  
  # ── 5. Return ─────────────────────────────────────────────────────────────────
  list(
    error_total = error_total,
    resid_time  = resid_time
  )
}

run_fpca_systems <- function(df,
                             group_var,
                             time_vars,
                             trait_col,
                             ti_col = NULL,
                             ...) {
  
  results <- lapply(names(time_vars), function(sys) {
    cat(sprintf("\nRunning FPCA for system: %s\n", sys))
    run_fpca_workflow(
      df        = df,
      group_var = group_var,
      time_var  = time_vars[[sys]],
      trait_col = trait_col,
      ti_col    = ti_col,
      ...
    )
  })
  names(results) <- names(time_vars)
  results
}

comparison <- run_fpca_systems(
  df        = my_dat1,
  group_var = "un_num",
  time_vars = c(abs_rnk = "Leaf", rel_ti = "rel_t", rel_e = "rel_e", rel_n = "rel_n"),
  trait_col = "Length"
)


extract_fpca_features <- function(comparison,
                                  sys = "rel_n",
                                  meta_df,
                                  group_var = "un_num",
                                  meta_vars = c("genotype", "origin",
                                                "LN_total", "LN_ear", "t_i"),
                                  min_time = 0.25) {
  
  res <- comparison[[sys]]
  fit <- res$fit
  
  # ── 1. Scores ───────────────────────────────────────────────────────────────
  scores <- as.data.frame(fit$xiEst) %>%
    setNames(paste0("fpc", seq_len(ncol(.)))) %>%
    dplyr::mutate(subject = res$subject_ids,
                  system  = sys)
  
  # ── 2. Metadata ─────────────────────────────────────────────────────────────
  meta <- meta_df %>%
    dplyr::select(all_of(c(group_var, meta_vars))) %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::distinct() %>%
    dplyr::rename(subject = all_of(group_var))
  
  # ── 3. Reconstructed curves + metadata ─────────────────────────────────────
  recon_meta <- res$recon_df %>%
    dplyr::left_join(meta, by = "subject") %>%
    dplyr::filter(time >= min_time)
  
  # ── 4. Longest leaf + canopy centre of mass per plant ─────────────────────
  longest_leaf <- recon_meta %>%
    dplyr::group_by(subject, across(all_of(meta_vars))) %>%
    dplyr::group_modify(~ {
      
      ln_total <- .y$LN_total
      
      # Convert rel_n -> leaf number, round to nearest integer,
      # average recon within each rounded leaf, then reconvert to rel_n
      by_leaf <- .x %>%
        dplyr::mutate(leaf_num_rounded = round(time * ln_total)) %>%
        dplyr::group_by(leaf_num_rounded) %>%
        dplyr::summarise(
          recon       = mean(recon, na.rm = TRUE),
          rel_n_whole = leaf_num_rounded / ln_total,
          .groups     = "drop"
        )
      
      # Longest leaf: max of discrete averaged positions
      best_row <- by_leaf[which.max(by_leaf$recon), ]
      
      # Canopy CoM: weighted average leaf position, weights = leaf length
      com_pos <- sum(by_leaf$rel_n_whole * by_leaf$recon) /
        sum(by_leaf$recon)
      
      tibble::tibble(
        longest_leaf     = best_row$recon,
        longest_leaf_pos = best_row$rel_n_whole,
        rate             = best_row$recon / best_row$rel_n_whole,
        com_pos          = com_pos
      )
    }) %>%
    dplyr::ungroup()
  
  # ── 5. Ascending limb derivatives ──────────────────────────────────────────
  derivs <- recon_meta %>%
    dplyr::group_by(subject, across(all_of(meta_vars))) %>%
    dplyr::group_modify(~ {
      
      peak_idx <- which.max(.x$recon)
      asc      <- .x[seq_len(peak_idx), ]
      
      dx <- diff(asc$time)
      dy <- diff(asc$recon)
      
      tibble::tibble(
        rank       = asc$time[-nrow(asc)],
        blade_size = asc$recon[-nrow(asc)],
        dy_dx      = dy / dx
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(env      = origin,
                  group_id = subject)
  
  # ── 6. Return ───────────────────────────────────────────────────────────────
  list(
    scores       = scores,
    recon        = recon_meta,
    longest_leaf = longest_leaf,
    derivs       = derivs
  )
}
# ==============================================================================
# 3. Extract values ----------------
# ==============================================================================


features <- extract_fpca_features(
  comparison = comparison,
  sys        = "rel_n",
  meta_df    = my_dat1,
  group_var  = "un_num"
)

# access outputs
features$scores
features$recon
features$longest_leaf
features$derivs       # feeds directly into find_max_slope_window



score_cor_table <- function(features, fpc_cols = NULL) {
  
  # ── 1. Join scores and longest leaf ─────────────────────────────────────────
  df <- features$scores %>%
    dplyr::left_join(features$longest_leaf, by = "subject")
  
  # ── 2. Identify FPC columns ──────────────────────────────────────────────────
  if (is.null(fpc_cols)) {
    fpc_cols <- grep("^fpc", names(df), value = TRUE)
  }
  
  # ── 3. Traits to correlate against ──────────────────────────────────────────
  trait_cols <- c("longest_leaf", "longest_leaf_pos", "com_pos")
  
  # ── 4. Compute correlations and p-values ────────────────────────────────────
  purrr::map_dfr(trait_cols, function(trait) {
    purrr::map_dfr(fpc_cols, function(fpc) {
      x <- df[[fpc]]
      y <- df[[trait]]
      complete <- complete.cases(x, y)
      
      if (sum(complete) < 3) {
        return(tibble::tibble(
          trait = trait, fpc = fpc,
          r = NA_real_, p = NA_real_, n = sum(complete)
        ))
      }
      
      ct <- cor.test(x[complete], y[complete], method = "pearson")
      tibble::tibble(
        trait = trait,
        fpc   = fpc,
        r     = round(ct$estimate, 3),
        p     = round(ct$p.value, 4),
        n     = sum(complete)
      )
    })
  }) %>%
    dplyr::mutate(
      sig = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ ""
      )
    )
}

cor_table <- score_cor_table(features)
cor_table


# #across multiple systems
# purrr::map_dfr(c("rel_n", "rel_ti", "rel_e", "abs_rnk"), function(sys) {
#   feat <- extract_fpca_features(comparison, sys = sys, meta_df = my_dat1)
#   score_cor_table(feat) %>% dplyr::mutate(system = sys)
# })



