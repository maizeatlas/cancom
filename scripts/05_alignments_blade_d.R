rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# Functional PCA of blade-length profiles: registration system comparison
# across platform and field experiments
#
#   Pipeline:
#     1. FPCA (fdapace) is run on per-plant blade-length curves under four
#        registration systems — absolute leaf rank, T.I.-relative rank,
#        ear-relative rank, and total-leaf-number normalised rank — for field (un_num) experiments.
#     2. Residual fit quality is benchmarked across registration systems via
#        total absolute residuals per plant and residual-by-time smooths,
#        enabling selection of the best-aligning time axis.
#     3. Covariance and correlation matrices are extracted per system;
#        off-diagonal correlations are summarised as mean correlation vs. lag
#        distance to characterise smoothness of the functional signal under
#        each registration.
#
#   Outputs:
#     A. CSV of field FPC scores + T.I. per genotype × environment
#        (05_dat_field_blade_fpc.csv).
#     B. Fig panel A — population mean curves per registration system
#        (blade length vs. registration coordinates), with T.I. range shaded.
#     C. Fig panel B — correlation heatmaps per registration system.
#     D. Fig panel C — variance profiles along the registration grid.
#     E. Fig panel D — mean off-diagonal correlation vs. lag distance.
#     F. Combined three-panel figure (A | B | C) saved as PDF.
# ==============================================================================


library(tidyverse) #data manipulation
library(data.table) #data manipulation
library(zoo)#data manipulation
library(ggh4x)#combining plots
library(patchwork) #combining plots
library(cowplot)#combining plots
library(ggridges)# plots
library(fdapace) #fda analysis
library(gtable)# Calculate approximate space needed:
library(grid)# Calculate approximate space needed:
library(ggnewscale)

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
run_fpca <- function(df, group_var, time_var, trait_col, ti_col) {
  
  # Filter NAs for both trait and time variable
  df <- df %>% 
    dplyr::filter(!is.na(.data[[trait_col]]), 
                  !is.na(.data[[time_var]]),
                  !is.na(.data[[ti_col]])) %>% 
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
  
  ti_vector <- df %>%
    group_by(.data[[group_var]]) %>%
    summarise(ti = first(.data[[ti_col]]), .groups = "drop") %>%
    pull(ti)
  
  attr(result, "group_ids") <- group_ids
  attr(result, "group_var") <- group_var
  attr(result, "t_i") <- ti_vector
  attr(result, "ti_col") <- ti_col
  
  return(result)
  
}


# Field data (blade)
fe_b_l <- run_fpca(my_dat1, group_var = "un_num", time_var = "Leaf", trait_col = "Length", ti_col = "t_i")   
fe_b_t <- run_fpca(my_dat1, group_var = "un_num", time_var = "rel_t", trait_col = "Length", ti_col = "t_i")
fe_b_e <- run_fpca(my_dat1, group_var = "un_num", time_var = "rel_e", trait_col = "Length", ti_col = "t_i")
fe_b_n <- run_fpca(my_dat1, group_var = "un_num", time_var = "rel_n", trait_col = "Length", ti_col = "t_i")

# ==============================================================================
# 3a. Extract results
# ==============================================================================

process_fpca_results <- function(fpca_obj, obj_name, trait_name,
                                 join_data, join_by, group_by_vars) {
  

  group_var <- attr(fpca_obj, "group_var")
  
  # base table: one row per plant curve
  base_tbl <- tibble(
    trait_name = trait_name,
    !!group_var := attr(fpca_obj, "group_ids"),
    fpca_source = obj_name,
    t_i = attr(fpca_obj, "t_i")
  )
  
  # FPCA scores
  scores_tbl <- as.data.frame(fpca_obj$xiEst) %>%
    setNames(paste0("fpc", seq_len(ncol(.))))
  
  base_tbl %>%
    bind_cols(scores_tbl) %>%
    left_join(join_data, by = setNames(join_by, group_var)) %>%
    
    # collapse plants -> genotype × origin
    group_by(across(all_of(c("trait_name", "fpca_source", group_by_vars)))) %>%
    
    summarise(
      across(starts_with("fpc"), ~ mean(.x, na.rm = TRUE)),
      t_i = mean(t_i, na.rm = TRUE),   # <- now a proper genotype×origin phenotype
      n_plants = n(),
      .groups = "drop"
    )
}


# Prepare join data
field_join_data <- my_dat1 %>% 
  dplyr::select(origin, un_num, range, row, block, genotype) %>% 
  distinct_all()

# Process all field FPCA results
field_results <- bind_rows(
  process_fpca_results(fe_b_l, "fe_b_l", "blade_length", field_join_data, "un_num", c("origin", "genotype")),
  process_fpca_results(fe_b_t, "fe_b_t", "blade_length", field_join_data, "un_num", c("origin", "genotype")),
  process_fpca_results(fe_b_e, "fe_b_e", "blade_length", field_join_data, "un_num", c("origin", "genotype")),
  process_fpca_results(fe_b_n, "fe_b_n", "blade_length", field_join_data, "un_num", c("origin", "genotype"))
)

# Convert to long format
field_results_long <- field_results %>%
  pivot_longer(
    cols = -c(trait_name, fpca_source, origin, genotype, n_plants, t_i),  # Pivot everything EXCEPT these
    names_to = "trait", 
    values_to = "value"
  )



# Function to extract variances (diagonal) and covariances (upper triangle)
extract_values <- function(mat, mat_name, coords = NULL) {
  # Convert to matrix if it's a list or data frame
  if (is.list(mat) && !is.matrix(mat)) {
    mat <- as.matrix(mat)
  } else if (is.data.frame(mat)) {
    mat <- as.matrix(mat)
  }
  # Ensure it's numeric
  mat <- matrix(as.numeric(mat), nrow = nrow(mat), ncol = ncol(mat))
  
  # Set default coordinates if not provided
  if (is.null(coords)) {
    coords <- 1:nrow(mat)
  }
  
  # Check that coords length matches matrix dimension
  if (length(coords) != nrow(mat)) {
    stop("Length of coords must equal the number of rows/columns in the matrix")
  }
  
  # Variances (diagonal)
  diag_indices <- 1:nrow(mat)
  variances <- data.frame(
    row_coord = coords[diag_indices],
    col_coord = coords[diag_indices],
    value = diag(mat),
    type = "Variance",
    matrix = mat_name
  )
  
  # upper triangle, excluding diagonal
  upper_indices <- which(upper.tri(mat), arr.ind = TRUE)
  covariances <- data.frame(
    row_coord = coords[upper_indices[, "row"]],
    col_coord = coords[upper_indices[, "col"]],
    value = mat[upper_indices],
    type = "Covariance",
    matrix = mat_name
  )
  
  rbind(variances, covariances)
}

# Helper function to process multiple matrices
process_matrices <- function(mat_list, coord_list, mat_names, type_label = NULL) {
  # Combine all matrices
  combined <- do.call(rbind, Map(extract_values, mat_list, mat_names, coord_list))
  
  # Convert matrix to factor with specific order
  combined$matrix <- factor(combined$matrix, levels = mat_names)
  
  # Handle type column
  if (!is.null(type_label)) {
    combined$type <- as.factor(type_label)
  } else {
    combined$type <- factor(combined$type, levels = c("Variance", "Covariance"))
  }
  
  return(combined)
}


# ==============================================================================
# 3b. Residual check by system----------------
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


# ==============================================================================
# 3d. Field Exp. -- Extract results
# ==============================================================================

# coords
fe_l_coords <- fe_b_l$workGrid
fe_t_coords <- fe_b_t$workGrid
fe_e_coords <- fe_b_e$workGrid
fe_n_coords <- fe_b_n$workGrid

# splines
fe_mu_l <- data.frame(coords = fe_l_coords, mu = fe_b_l$mu)  
fe_mu_t <- data.frame(coords = fe_t_coords, mu = fe_b_t$mu)
fe_mu_e <- data.frame(coords = fe_e_coords, mu = fe_b_e$mu)  
fe_mu_n <- data.frame(coords = fe_n_coords, mu = fe_b_n$mu)

# covariance
fe_mat_l_cov <- fe_b_l$fittedCov
fe_mat_t_cov <- fe_b_t$fittedCov       
fe_mat_e_cov <- fe_b_e$fittedCov   
fe_mat_n_cov <- fe_b_n$fittedCov  

# correlation
fe_mat_l_cor <- cov2cor(as.matrix(fe_mat_l_cov))
fe_mat_t_cor <- cov2cor(as.matrix(fe_mat_t_cov))
fe_mat_e_cor <- cov2cor(as.matrix(fe_mat_e_cov))
fe_mat_n_cor <- cov2cor(as.matrix(fe_mat_n_cov))

# Process covariance matrices
fe_combined_cov <- process_matrices(
  mat_list = list(fe_mat_l_cov, fe_mat_t_cov, fe_mat_e_cov, fe_mat_n_cov),
  coord_list = list(fe_l_coords, fe_t_coords, fe_e_coords, fe_n_coords),
  mat_names = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm."),
)

# Process correlation matrices
fe_combined_cor <- process_matrices(
  mat_list = list(fe_mat_l_cor, fe_mat_t_cor, fe_mat_e_cor, fe_mat_n_cor),
  coord_list = list(fe_l_coords, fe_t_coords, fe_e_coords, fe_n_coords),
  mat_names = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm."),
  type_label = NULL  # Don't override the type, keep "Variance" and "Covariance"
) %>%
  mutate(type = case_when(
    type == "Variance" ~ "Cor_diag",
    type == "Covariance" ~ "Cor_offdiag",
    TRUE ~ type
  )) %>%
  group_by(matrix) %>%
  mutate(
    # Create index position for each coordinate
    row_idx = match(row_coord, sort(unique(c(row_coord, col_coord)))),
    col_idx = match(col_coord, sort(unique(c(row_coord, col_coord))))
  ) %>%
  ungroup()

# call residual##################
fe_comparison <- compare_fpca_systems(
  df        = my_dat1,
  group_var = "un_num",
  time_vars = list(abs_rnk = "Leaf", rel_ti = "rel_t", rel_e = "rel_e", rel_n = "rel_n"),
  trait_col = "Length"
)
# access directly — no intermediate steps needed
fe_comparison$error_total  # feeds p_error_box

resid_time<-fe_comparison$resid_time%>%
  dplyr::rename(matrix=system)
# feeds p_resid_smooth

# ==============================================================================
# 4. Plots 
# ==============================================================================

# Define "fe" ranges if TI for each matrix (CUSTOMIZE THESE!)
fe_leaf_sequence <- my_dat1 %>%
  group_by(un_num) %>%
  reframe(Leaf = 1:max(Leaf, na.rm = TRUE))

fe_rect_data <- fe_leaf_sequence %>%
  left_join(my_dat1, by = c("un_num", "Leaf")) %>%
  group_by(un_num) %>%
  fill(LN_total, t_i, .direction = "downup") %>%
  filter(!is.na(LN_total)) %>%
  mutate(
    ti_floor = as.integer(floor(t_i)),
    rel_t = -(ti_floor - Leaf),
    first_val = first(na.omit(rel_t)),
    first_pos = which(!is.na(rel_t))[1],
    rel_t = first_val + (row_number() - first_pos),
    rel_e = na.approx(rel_e, na.rm = FALSE, rule = 2),
    rel_n = na.approx(rel_n, na.rm = FALSE, rule = 2),
    maxlr = max(t_i, na.rm = TRUE),
    minlr = min(t_i, na.rm = TRUE)
  ) %>%
  slice_min(abs(rel_t), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  summarize(
    t_min_l = min(minlr, na.rm = TRUE),
    t_max_l = max(maxlr, na.rm = TRUE),
    t_min_t = -0.5, # min(rel_t, na.rm = TRUE),
    t_max_t = 0.5,  # max(rel_t, na.rm = TRUE),
    t_min_e = min(rel_e, na.rm = TRUE),
    t_max_e = max(rel_e, na.rm = TRUE),
    t_min_n = min(rel_n, na.rm = TRUE),
    t_max_n = max(rel_n, na.rm = TRUE)
  ) %>%
  {data.frame(
    matrix = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm."),
    xmin = c(.$t_min_l, .$t_min_t, .$t_min_e, .$t_min_n),
    xmax = c(.$t_max_l, .$t_max_t, .$t_max_e, .$t_max_n)
  )} %>%
  mutate(
    ymin = xmin,
    ymax = xmax,
    matrix = factor(matrix, levels = matrix)
  )

# Step 1: build per-un_num values, long format instead of summarized
fe_density_long <- fe_leaf_sequence %>%
  left_join(my_dat1, by = c("un_num", "Leaf")) %>%
  group_by(un_num) %>%
  fill(LN_total, t_i, .direction = "downup") %>%
  filter(!is.na(LN_total)) %>%
  mutate(
    ti_floor = as.integer(floor(t_i)),
    rel_t = -(ti_floor - Leaf),
    first_val = first(na.omit(rel_t)),
    first_pos = which(!is.na(rel_t))[1],
    rel_t = first_val + (row_number() - first_pos),
    rel_e = na.approx(rel_e, na.rm = FALSE, rule = 2),
    rel_n = na.approx(rel_n, na.rm = FALSE, rule = 2),
    maxlr = max(t_i, na.rm = TRUE),
    minlr = min(t_i, na.rm = TRUE)
  ) %>%
  slice_min(abs(rel_t), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  pivot_longer(
    cols = c(minlr, maxlr, rel_t, rel_e, rel_n),
    names_to = "var_key",
    values_to = "value"
  ) %>%
  mutate(
    matrix = case_when(
      var_key %in% c("minlr", "maxlr") ~ "Abs. Rank",
      var_key == "rel_t"               ~ "Rel. T.I.",
      var_key == "rel_e"               ~ "Rel. Ear Rank",
      var_key == "rel_n"               ~ "Leaf No. Norm."
    ),
    matrix = factor(matrix, levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm."))
  ) %>%
  filter(is.finite(value))

# Step 2: density per matrix category
rect_grad_fe <- fe_density_long %>%
  group_by(matrix) %>%
  group_modify(~ {
    if (nrow(.x) < 2) {
      return(data.frame(xmin = numeric(0), xmax = numeric(0),
                        fill = numeric(0), x = numeric(0), dens = numeric(0)))
    }
    d <- density(.x$value)
    data.frame(
      xmin = head(d$x, -1),
      xmax = tail(d$x, -1),
      fill = head(d$y, -1),
      x    = head(d$x, -1),
      dens = head(d$y, -1)
    )
  }) %>%
  group_by(matrix) %>%
  mutate(
    fill_norm = fill / max(fill, na.rm = TRUE),
    dens_norm = dens / max(dens, na.rm = TRUE)
  ) %>%
  ungroup()

# Step 3: scale dens line to a fraction of each panel's y-range
ylims_fe <- c(`Abs. Rank` = 85, `Rel. T.I.` = 85, `Rel. Ear Rank` = 85, `Leaf No. Norm.` = 85)  # fill in real limits

rect_grad_fe <- rect_grad_fe %>%
  mutate(dens_scaled = dens_norm * 0.15 * ylims_fe[as.character(matrix)])


# ==============================================================================
# 4a. Field Exp. (fe) - Mean curves
# ==============================================================================
# Combine all data frames with a grouping variable
fe_mu_combined <- bind_rows(
  fe_mu_l %>% mutate(matrix = "Abs. Rank"),
  fe_mu_t %>% mutate(matrix = "Rel. T.I."),
  fe_mu_e %>% mutate(matrix = "Rel. Ear Rank"),
  fe_mu_n %>% mutate(matrix = "Leaf No. Norm.")
)

fe_mu_combined <- fe_mu_combined %>%
  mutate(matrix = factor(matrix, levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm.")))

# Create custom x-scales for each facet
x_scales_mu <- fe_mu_combined %>%
  distinct(matrix, coords) %>%
  arrange(matrix, coords) %>%
  split(.$matrix) %>%
  lapply(function(df) {
    
    m <- unique(df$matrix)
    
    coord_limits <- c(min(df$coords), max(df$coords))
    
    if (m %in% c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank")) {
      # Grid lines at whole numbers
      coord_grid <- seq(ceiling(coord_limits[1]), floor(coord_limits[2]), by = 1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 0)
      
    } else if (m == "Leaf No. Norm.") {
      # Grid lines at 0.1 intervals
      coord_grid <- seq(ceiling(coord_limits[1] * 10) / 10, 
                        floor(coord_limits[2] * 10) / 10, 
                        by = 0.1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
      
    } else {
      coord_grid <- pretty(coord_limits, n = 5)
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
    }
    
    scale_x_continuous(
      breaks = coord_label_at,      # Only show labels at min/max
      labels = coord_labels,
      limits = coord_limits,
      minor_breaks = coord_grid,    # Grid lines at whole numbers or 0.1
      expand = expansion(mult = c(0.1, 0.1))
    )
  })

# Create faceted plot
p_fe_mu <- ggplot(fe_mu_combined, aes(x = coords, y = mu, color = matrix)) +
  # shading before lines
  geom_rect(data = rect_grad_fe,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_norm),
            inherit.aes = FALSE) +
  scale_fill_gradient(low = alpha("grey85", 0.2), high = alpha("grey50", 0.7), guide = "none") +
  facet_wrap(~ matrix, scales = "free") +
  geom_line(linewidth = 1.2) +
  ggh4x::facet_wrap2(~ matrix, scales = "free_x", ncol = 4, axes = "all") +
  ggh4x::facetted_pos_scales(
    x = x_scales_mu
  ) +
  scale_color_manual(values = c(
    "Abs. Rank"      = "#0072B2",
    "Rel. Ear Rank"  = "#D55E00",
    "Rel. T.I."      = "#009E73",
    "Leaf No. Norm." = "#CC79A7"
  )) +
  labs(title = "E", x = "Registration coordinates", y = "Blade length\n(cm)") +
  scale_y_continuous(limits = c(20, 105),
                     breaks = function(lim) seq(0, 100, 20)) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),  # Turn off major grid
    panel.grid.minor.x = element_line(linewidth = 0.3, colour = "grey90"),  # Use minor for custom grid
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(),
    panel.spacing = unit(1, "lines")
  )
p_fe_mu

# ==============================================================================
# 4b. Field Exp. (fe) - Correlation matrices (heat map)
# ==============================================================================
# Convert long format to heatmap format
# For upper triangle + diagonal, we need to create the full matrix visualization
fe_heatmap_data <- fe_combined_cor %>%
  dplyr::mutate(matrix = factor(matrix, levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm.")))

# Rename columns for the plot
fe_heatmap_data <- fe_heatmap_data %>%
  dplyr::rename(row = row_coord, col = col_coord)

# Create custom x-scales for each facet
x_scales_cor <- fe_heatmap_data %>%
  distinct(matrix, col) %>%
  dplyr::arrange(matrix, col) %>%
  split(.$matrix) %>%
  lapply(function(df) {
    
    m <- unique(df$matrix)
    
    coord_limits <- c(min(df$col), max(df$col))
    
    if (m %in% c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank")) {
      # Grid lines at whole numbers
      coord_grid <- seq(ceiling(coord_limits[1]), floor(coord_limits[2]), by = 1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 0)
      
    } else if (m == "Leaf No. Norm.") {
      # Grid lines at 0.1 intervals
      coord_grid <- seq(ceiling(coord_limits[1] * 10) / 10, 
                        floor(coord_limits[2] * 10) / 10, 
                        by = 0.1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
      
    } else {
      coord_grid <- pretty(coord_limits, n = 5)
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
    }
    
    scale_x_continuous(
      breaks = coord_label_at,      # Only show labels at min/max
      labels = coord_labels,
      limits = coord_limits,
      minor_breaks = coord_grid,    # Grid lines at whole numbers or 0.1
      expand = expansion(mult = c(0.1, 0.1))
    )
  })

# Create custom y-scales for each facet
y_scales_cor <- fe_heatmap_data %>%
  distinct(matrix, row) %>%
  dplyr::arrange(matrix, row) %>%
  split(.$matrix) %>%
  lapply(function(df) {
    
    m <- unique(df$matrix)
    
    coord_limits <- c(min(df$row), max(df$row))
    
    if (m %in% c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank")) {
      # Grid lines at whole numbers
      coord_grid <- seq(ceiling(coord_limits[1]), floor(coord_limits[2]), by = 1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 0)
      
    } else if (m == "Leaf No. Norm.") {
      # Grid lines at 0.1 intervals
      coord_grid <- seq(ceiling(coord_limits[1] * 10) / 10, 
                        floor(coord_limits[2] * 10) / 10, 
                        by = 0.1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
      
    } else {
      coord_grid <- pretty(coord_limits, n = 5)
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
    }
    
    scale_y_reverse(
      breaks = coord_label_at,      # Only show labels at min/max
      labels = coord_labels,
      limits = rev(coord_limits),   # Reversed for scale_y_reverse
      minor_breaks = coord_grid,    # Grid lines at whole numbers or 0.1
      expand = expansion(mult = c(0.1, 0.1))
    )
  })

# Create heatmap with highlighted regions and matching colors
p_fe_cor <- ggplot(fe_heatmap_data, aes(x = col, y = row, fill = value)) +
  geom_raster() +
  geom_rect(data = fe_rect_data,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = NA, color = "grey", linewidth = 1, inherit.aes = FALSE) +
  ggh4x::facet_wrap2(~ matrix, ncol = 4, scales = "free", axes = "all") +
  ggh4x::facetted_pos_scales(
    x = x_scales_cor,
    y = y_scales_cor
  ) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, 
    limits = c(-1, 1),                 # optional but recommended
    breaks = seq(-1, 1, by = 0.5),      # <-- 0.5 breaks
    name = "Correlation",
    guide = guide_colorbar(
      barwidth = unit(3, "cm"),
      barheight = unit(0.5, "cm"),
      title.position = "left",
      title.vjust = 0.8
    )
  ) +
  labs(title = "F", x = "Registration coordinates", y = "Registration\ncoordinates") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    # legend.position = c(0.65, 1.5),
    # legend.justification = c(0, 1),
    # legend.direction = "horizontal",
    # legend.box.spacing = unit(0, "cm"),
    # legend.text = element_text(size = 12),
    # legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    # legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),  # Turn off major grid
    panel.grid.minor.x = element_line(linewidth = 0.3, colour = "grey90"),  # Use minor for custom grid
    panel.grid.major.y = element_blank(),  # Turn off major grid
    panel.grid.minor.y = element_line(linewidth = 0.3, colour = "grey90"),  # Use minor for custom grid
    strip.text = element_text(),
    panel.spacing = unit(1, "lines")
  )
p_fe_cor

# # Version without legend
p_fe_cor_no_legend <- p_fe_cor +
  theme(legend.position = "none")
p_fe_cor_no_legend

# ==============================================================================
# 4c. Field Exp. (fe) - Scaled variance across grid (diagonal)
# ==============================================================================
# Plot variance by rank position
fe_combined_var <- fe_combined_cov %>%
  dplyr::filter(type == "Variance") %>%
  dplyr::group_by(matrix) %>%
  dplyr::mutate(idx = row_number(),
                value_s = as.vector(scale(value, center=T, scale=T)))

fe_combined_var <- fe_combined_var %>%
  mutate(matrix = factor(matrix, levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm.")))

# Create custom x-scales for each facet
x_scales_var <- fe_combined_var %>%
  distinct(matrix, row_coord) %>%
  arrange(matrix, row_coord) %>%
  split(.$matrix) %>%
  lapply(function(df) {
    
    m <- unique(df$matrix)
    
    coord_limits <- c(min(df$row_coord), max(df$row_coord))
    
    if (m %in% c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank")) {
      # Grid lines at whole numbers
      coord_grid <- seq(ceiling(coord_limits[1]), floor(coord_limits[2]), by = 1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 0)
      
    } else if (m == "Leaf No. Norm.") {
      # Grid lines at 0.1 intervals
      coord_grid <- seq(ceiling(coord_limits[1] * 10) / 10, 
                        floor(coord_limits[2] * 10) / 10, 
                        by = 0.1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
      
    } else {
      coord_grid <- pretty(coord_limits, n = 5)
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
    }
    
    scale_x_continuous(
      breaks = coord_label_at,      # Only show labels at min/max
      labels = coord_labels,
      limits = coord_limits,
      minor_breaks = coord_grid,    # Grid lines at whole numbers or 0.1
      expand = expansion(mult = c(0.1, 0.1))
    )
  })

# Create faceted plot
p_fe_var <- ggplot(fe_combined_var, aes(x = row_coord, y = value, color = matrix))  +
  # shading before lines
  # geom_rect(data = fe_rect_data,
  #           aes(xmin = xmin, xmax = xmax,
  #               ymin = -Inf, ymax = Inf),
  #           inherit.aes = FALSE,
  #           alpha = 1, fill = "grey") +
  geom_rect(data = rect_grad_fe,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_norm),
            inherit.aes = FALSE) +
  scale_fill_gradient(low = alpha("grey85", 0.2), high = alpha("grey50", 0.7), guide = "none") +
  facet_wrap(~ matrix, scales = "free") +
  geom_line(linewidth = 1.2) +
  ggh4x::facet_wrap2(~ matrix, scales = "free_x", ncol = 4, axes = "all") +
  ggh4x::facetted_pos_scales(
    x = x_scales_var
  ) +
  scale_color_manual(values = c(
    "Abs. Rank"      = "#0072B2",
    "Rel. Ear Rank"  = "#D55E00",
    "Rel. T.I."      = "#009E73",
    "Leaf No. Norm." = "#CC79A7"
  )) +
  labs(title = "G", x = "Registration coordinates", y = "Variance\n") +
  # scale_y_continuous(labels = scales::label_number(accuracy = 0.1),
  #                    limits = c(-2.0, 3.0),
  #                    breaks = function(lim) seq(-2.0, 3.0, 1)) +
  scale_y_continuous(limits = c(0, 325),
                     breaks = function(lim) seq(0, 350, 100)) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),  # Turn off major grid
    panel.grid.minor.x = element_line(linewidth = 0.3, colour = "grey90"),  # Use minor for custom grid
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(),
    panel.spacing = unit(1, "lines")
  )
p_fe_var

## summary for ranks < 4 phytomers from max (idx < 45,  ~26 ranks) 
fe_combined_var %>%
  filter(idx < 45) %>%
  group_by(matrix) %>%
  summarize(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  ) %>%
  mutate(
    pct_reduction = (max_value[matrix == "Abs. Rank"] - max_value) / max_value[matrix == "Abs. Rank"] * 100,
    fold_increase = max_value[matrix == "Abs. Rank"] / max_value
  )

## summary for ranks < 4 phytomers from max (idx < 39,  ~23 ranks) 
fe_combined_var %>%
  filter(idx < 39) %>%
  group_by(matrix) %>%
  summarize(
    min_value = min(value, na.rm = TRUE),
    max_value = max(value, na.rm = TRUE)
  ) %>%
  mutate(
    pct_reduction = (max_value[matrix == "Abs. Rank"] - max_value) / max_value[matrix == "Abs. Rank"] * 100,
    fold_increase = max_value[matrix == "Abs. Rank"] / max_value
  )

# # ==============================================================================
# # 4d. Field Exp. (fe) - Mean covariance across grid (off-diagonal)
# # ==============================================================================
# # Calculate mean correlation per lag distance
fe_combined_cor_lag <- fe_combined_cor %>%
  dplyr::filter(type == "Cor_offdiag") %>%
  mutate(lag_idx = abs(col_idx - row_idx),
         lag_rank = round(abs(col_coord - row_coord), 2)) %>%
  dplyr::group_by(matrix, lag_idx) %>%
  dplyr::summarise(mean_cov = mean(value, na.rm = TRUE),
            lag_rank = first(lag_rank),  # Or median(lag_rank), min(lag_rank), etc.
            .groups = "drop")

fe_combined_cor_lag <- fe_combined_cor_lag %>%
  mutate(matrix = factor(matrix, levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm.")))

rank_to_idx <- function(target_rank, lag_rank, lag_idx) {
  approx(
    x = lag_rank,
    y = lag_idx,
    xout = target_rank,
    rule = 2
  )$y
}

x_scales <- fe_combined_cor_lag %>%
  distinct(matrix, lag_idx, lag_rank) %>%
  dplyr::arrange(matrix, lag_idx) %>%
  split(.$matrix) %>%
  lapply(function(df) {
    
    m <- unique(df$matrix)
    
    rank_limits <- c(min(df$lag_rank), max(df$lag_rank))
    
    if (m %in% c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank")) {
      # Grid lines at whole numbers
      rank_grid <- seq(ceiling(rank_limits[1]), floor(rank_limits[2]), by = 1)
      # Labels only at min and max
      rank_label_at <- rank_limits
      rank_labels <- round(rank_limits, 0)
      
    } else if (m == "Leaf No. Norm.") {
      # Grid lines at 0.1 intervals
      rank_grid <- seq(ceiling(rank_limits[1] * 10) / 10, 
                       floor(rank_limits[2] * 10) / 10, 
                       by = 0.1)
      # Labels only at min and max
      rank_label_at <- rank_limits
      rank_labels <- round(rank_limits, 1)
      
    } else {
      rank_grid <- pretty(rank_limits, n = 5)
      rank_label_at <- rank_limits
      rank_labels <- round(rank_limits, 1)
    }
    
    # Convert grid positions to idx
    idx_grid <- rank_to_idx(
      target_rank = rank_grid,
      lag_rank = df$lag_rank,
      lag_idx  = df$lag_idx
    )
    
    # Convert label positions to idx
    idx_labels <- rank_to_idx(
      target_rank = rank_label_at,
      lag_rank = df$lag_rank,
      lag_idx  = df$lag_idx
    )
    
    idx_limits <- rank_to_idx(
      target_rank = rank_limits,
      lag_rank = df$lag_rank,
      lag_idx  = df$lag_idx
    )
    
    scale_x_continuous(
      breaks = idx_labels,      # Only show labels at min/max
      labels = rank_labels,
      limits = idx_limits,
      minor_breaks = idx_grid,  # Grid lines at whole numbers or 0.1
      expand = expansion(mult = c(0.1, 0.1))
    )
  })

# Create lag plot
p_fe_cov_lag <- ggplot(fe_combined_cor_lag, aes(x = lag_idx, y = mean_cov, color = matrix)) +
  geom_line(linewidth = 1.2) +
  scale_y_continuous(limits = c(-0.5, 1.0),
                     breaks = function(lim) seq(-0.5, 1.0, 0.5)) +
  ggh4x::facet_wrap2(~ matrix, scales = "free_x", ncol = 4, axes = "all") +
  ggh4x::facetted_pos_scales(
    x = x_scales
  ) +
  scale_color_manual(values = c(
    "Abs. Rank"      = "#0072B2",
    "Rel. Ear Rank"  = "#D55E00",
    "Rel. T.I."      = "#009E73",
    "Leaf No. Norm." = "#CC79A7"
  )) +
  labs(title = "H", x = "Registration lag distance", y = "Mean\ncorrelation") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),  # Turn off major grid
    panel.grid.minor.x = element_line(linewidth = 0.3, colour = "grey90"),  # Use minor for custom grid
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(),
    panel.spacing = unit(1, "lines")
  )
p_fe_cov_lag




# ==============================================================================
# 3. fe -residual
# ==============================================================================

# Create custom x-scales for each facet
x_scales_var <- fe_combined_var %>%
  distinct(matrix, row_coord) %>%
  arrange(matrix, row_coord) %>%
  split(.$matrix) %>%
  lapply(function(df) {
    
    m <- unique(df$matrix)
    
    coord_limits <- c(min(df$row_coord), max(df$row_coord))
    
    if (m %in% c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank")) {
      # Grid lines at whole numbers
      coord_grid <- seq(ceiling(coord_limits[1]), floor(coord_limits[2]), by = 1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 0)
      
    } else if (m == "Leaf No. Norm.") {
      # Grid lines at 0.1 intervals
      coord_grid <- seq(ceiling(coord_limits[1] * 10) / 10, 
                        floor(coord_limits[2] * 10) / 10, 
                        by = 0.1)
      # Labels only at min and max
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
      
    } else {
      coord_grid <- pretty(coord_limits, n = 5)
      coord_label_at <- coord_limits
      coord_labels <- round(coord_limits, 1)
    }
    
    scale_x_continuous(
      breaks = coord_label_at,      # Only show labels at min/max
      labels = coord_labels,
      limits = coord_limits,
      minor_breaks = coord_grid,    # Grid lines at whole numbers or 0.1
      expand = expansion(mult = c(0.1, 0.1))
    )
  })



p_resid_smooth <- ggplot(
  resid_time %>%
    dplyr::mutate(matrix = dplyr::recode(matrix,
                                         "abs_rnk" = "Abs. Rank",
                                         "rel_ti"  = "Rel. T.I.",
                                         "rel_e"   = "Rel. Ear Rank",
                                         "rel_n"   = "Leaf No. Norm."
    )) %>%
    dplyr::mutate(matrix = factor(matrix,
                                  levels = c("Abs. Rank", "Rel. T.I.", "Rel. Ear Rank", "Leaf No. Norm."))),
  aes(x = time, y = resid, colour = matrix)   # <-- fill = matrix removed from here
) +
  geom_rect(data = rect_grad_fe,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill_norm),
            inherit.aes = FALSE) +
  scale_fill_gradient(low = alpha("grey85", 0.2), high = alpha("grey50", 0.7), guide = "none") +
  
  facet_wrap(~ matrix, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "grey40") +
  geom_point(alpha = 0.05, colour = "black") +
  geom_smooth(aes(fill = after_scale(colour)),   # <-- ribbon fill derived from colour, not a new scale
              method = "loess", se = TRUE, alpha = 0.5, linewidth = 1.2) +
  ggh4x::facet_wrap2(~ matrix, scales = "free_x", ncol = 4, axes = "all") +
  ggh4x::facetted_pos_scales(x = x_scales_var) +
  
  scale_colour_manual(values = c(
    "Abs. Rank"      = "#0072B2",
    "Rel. Ear Rank"  = "#D55E00",
    "Rel. T.I."      = "#009E73",
    "Leaf No. Norm." = "#CC79A7"
  )) +
  labs(
    title = "",
    x     = "Registration coordinates",
    y     = "Residual\n(cm)"
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey90"),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(),
    panel.spacing = unit(1, "lines")
  )

p_resid_smooth


# ==============================================================================
# 4. fe - COMBINE PLOTS
# ==============================================================================

# Create the adjusted plots
fixed_left_margin <- 25
p_fe_mu_adj <- p_fe_mu +
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, fixed_left_margin, "pt"),
    axis.title.y = element_text(margin = margin(r = 8))
  )
p_fe_cor_no_legend_adj <- p_fe_cor_no_legend +
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, fixed_left_margin, "pt"),
    axis.title.y = element_text(margin = margin(r = 15))
  )
p_fe_var_adj <- p_fe_var +
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, fixed_left_margin, "pt"),
    axis.title.y = element_text(margin = margin(r = 5))
  )
p_fe_cov_lag_adj <- p_fe_cov_lag +
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, fixed_left_margin, "pt"),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15))
  )
p_resid_smooth_adj<-p_resid_smooth+
  theme(
    plot.margin = margin(5.5, 5.5, 5.5, fixed_left_margin, "pt"),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.title.x = element_text(margin = margin(t = 15))
  )



# Change titles to A, B, C, D
p_fe_mu_adj <- p_fe_mu_adj + labs(title = "A")
p_fe_cor_no_legend_adj <- p_fe_cor_no_legend_adj + labs(title = "B")
p_resid_smooth_adj  <-  p_resid_smooth_adj+ labs(title = "C")
p_fe_var_adj <- p_fe_var_adj + labs(title = "C")
p_fe_cov_lag_adj <- p_fe_cov_lag_adj + labs(title = "D")

# Convert to grobs
g_mu <- ggplotGrob(p_fe_mu_adj)
g_cor <- ggplotGrob(p_fe_cor_no_legend_adj)
g_resid <- ggplotGrob(p_resid_smooth_adj)
g_var <- ggplotGrob(p_fe_var_adj)
g_lag <- ggplotGrob(p_fe_cov_lag_adj)


# Find the panel columns (where the actual plot panels are)
panel_cols <- g_mu$layout$l[grep("panel", g_mu$layout$name)]

# Get the maximum width for each column position across all plots
all_grobs <- list(g_mu, g_cor, g_resid)
max_widths <- do.call(unit.pmax, lapply(all_grobs, function(g) g$widths))

# Apply the maximum widths to all grobs
for(i in seq_along(all_grobs)) {
  all_grobs[[i]]$widths <- max_widths
}

# Combine using patchwork with the aligned grobs
fe_combined_plot <- wrap_plots(
  all_grobs[[1]], 
  all_grobs[[2]], 
  all_grobs[[3]], 
  # all_grobs[[4]],
  ncol = 1,
  heights = c(1, 1, 1)
)


#save it
ggsave(
  filename = paste0(fig_dir, "Fig4_fe_standalone_blade_", timestamp, ".pdf"), 
  plot = fe_combined_plot, 
  width = 7.5,
  height = 7.5,
  dpi = 300
)


# Simple plot just for the legend
fe_legend_temp <- ggplot(fe_heatmap_data,
                         aes(x = col, y = row, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),                 # optional but recommended
    breaks = seq(-1, 1, by = 0.5),      # <-- 0.5 breaks
    name = "Correlation",
    guide = guide_colorbar(
      barwidth = unit(3.0, "cm"),
      barheight = unit(0.5, "cm"),
      title.position = "left",
      title.vjust = 0.8
    )
  ) +
  theme(legend.direction = "horizontal")

# Extract the legend
fe_cor_legend_h <- get_legend(fe_legend_temp)
fe_cor_legend_h <- ggdraw(fe_cor_legend_h)

# # View it
fe_cor_legend_h


#save it
ggsave(
  filename = paste0(fig_dir, "fe_cor_legend_h_", timestamp, ".pdf"), 
  plot = fe_cor_legend_h, 
  width = 2.25, 
  height = 0.5,  # Increased height for 4 stacked plots
  dpi = 300
)


