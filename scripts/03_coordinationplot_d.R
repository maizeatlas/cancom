rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# Phytomer elongation dynamics: blade, sheath, and internode growth relative
# to tassel initiation (T.I.)
#
#   Pipeline:
#     1. Linear interpolation of thermal time (tt20_bp) against visible-leaf
#        (phyllo) and ligule (colla) counts per plant yields phytomer-level
#        thermal intersection points across a continuous leaf grid (0.1 steps).
#     2. T.I. thermal time is extracted per plant by projecting the visible-leaf
#        regression onto the observed T.I. value.
#     3. P-spline curves (GAM, cubic spline basis, k = 11) are fitted per plant
#        for blade length, sheath length, and internode distance across the
#        phytomer axis; predictions are min-max normalized within plant.
#     4. Flexposure (thermal time relative to T.I.) is computed by subtracting
#        the T.I. thermal intercept from each phytomer intersection point, then
#        merged with the P-spline predictions matching phyllo/colla parameter
#        to trait type (blade → phyllo; sheath, internode → colla).
#
#   Outputs:
#     A. Fig3A — absolute length (cm) vs. flexposure, faceted by phytomer
#        component (fixed and free-y versions), with per-plant trajectories
#        and population GAM smooths.
#     B. Fig3B — normalized length vs. flexposure, overlaying all three
#        components on a common [0, 1] scale to reveal relative expansion
#        timing across blade, sheath, and internode.
# ==============================================================================

library(tidyverse) #data manipulation and ggplot
library(mgcv) #GAM
library(ggh4x) #facetting


# ---- Paths ----
dat_dir = paste0(getwd(),"/data/")
fig_dir = paste0(getwd(),"/figures/")

# ==============================================================================
# 1. Read and summarize data 
# ==============================================================================

#canopy dimensions
my_dat0 <- read.csv(paste0(dat_dir,"04_dat_platform_canopy.csv")) %>%
  dplyr::group_by(pot)%>%
  dplyr::mutate(t_i=(tot_leaf-6.5)/1.51)%>%
  dplyr::mutate(pot=as.factor(pot))%>%
  dplyr::filter(genotype!="CML258")%>%
  ungroup()


#phenological data
my_dat1 <- read.csv(paste0(dat_dir,"05_dat_platform_allpheno.csv")) %>%
  dplyr::mutate(pot=as.factor(pot))%>%
  left_join(my_dat0%>%
              dplyr::select(t_i,earleaf, pot)%>%
              distinct(), by="pot")%>%
  dplyr::filter(genotype!="CML258")

# ==============================================================================
# 2. Predict leaf tip and collar appearance based on platform observations
# ==============================================================================

compute_leaf_intersections <- function(data, x_var, y_var, 
                                       n_leaves = 22, 
                                       increment = 0.1,
                                       subtract_from_x = 0) {
  
  horizontal_lines <- seq(1, n_leaves, by = increment)
  
  data %>%
    dplyr::group_by(pot) %>%
    dplyr::reframe(
      intersection_point = {
        # Fit tt20_bp as a function of (visi - subtract_from_y)
        x_adjusted <- .data[[x_var]] - subtract_from_x
        fit <- lm(.data[[y_var]] ~ x_adjusted)
        
        # Predict tt20_bp for each leaf number
        predict(fit, newdata = data.frame(x_adjusted = horizontal_lines))
      },
      leaf = horizontal_lines
    ) %>%
    dplyr::arrange(pot, intersection_point) %>%
    dplyr::mutate(wleaf = floor(leaf)) %>%
    dplyr::left_join(
      my_dat1 %>% dplyr::select(pot, genotype, treatment) %>% distinct(),
      by = "pot"
    ) %>%
    dplyr::group_by(genotype, treatment, if(increment == 1) leaf else NULL) %>%
    dplyr::mutate(rep = as.factor(row_number())) %>%
    dplyr::ungroup() %>%
    droplevels()
}

# Phyllo: continuous leaves (0.1 increment)
int_phyllo <- compute_leaf_intersections(
  my_dat1, 
  x_var = "visi",      # Predictor
  y_var = "tt20_bp",   # Response
  n_leaves = 22,
  increment = 0.1,
  subtract_from_x = 0
) %>%
  dplyr::mutate(parameter="phyllo")

# Collar: discrete leaves (1 increment)
int_colla <- compute_leaf_intersections(
  my_dat1,
  x_var = "ligul",     # Predictor
  y_var = "tt20_bp",   # Response
  n_leaves = 22,
  increment = 0.1,
  subtract_from_x = 0
)%>%
  dplyr::mutate(parameter="colla")

# ==============================================================================
# 3. T.I. in thermal time
# ==============================================================================

# Separate thermal intersection function
compute_thermal_intersection <- function(data, x_var, y_var, thermal_var) {
  
  data %>%
    dplyr::group_by(pot) %>%
    dplyr::summarise(
      thrmal_intrsect = {
        # Extract vectors explicitly
        x_vals <- .data[[x_var]]
        y_vals <- .data[[y_var]]
        
        # Fit with explicit vectors
        fit <- lm(y_vals ~ x_vals)
        
        thermal_val <- unique(.data[[thermal_var]])
        if(length(thermal_val) > 1) thermal_val <- thermal_val[1]
        
        # Predict using coefficient names from the model
        new_df <- data.frame(x_vals = thermal_val)
        as.numeric(predict(fit, newdata = new_df))
      },
      .groups = "drop"
    ) %>%
    distinct() %>%
    dplyr::left_join(
      data %>% 
        dplyr::select(pot, genotype) %>% 
        distinct(),
      by = "pot"
    ) %>%
    droplevels()
}

# Usage
thermal_ti <- compute_thermal_intersection(my_dat1, 
   x_var = "visi",      
   y_var = "tt20_bp",   
   thermal_var = "t_i"  
)

# ==============================================================================
# 4. P-spline phytomer curves
# ==============================================================================

create_multiple_pspline <- function(
    data,
    y_vars,
    grid,
    k = 11) 
  {
  
  create_single_pspline <- function(data, y_var, grid, k) {
    # Create the fitted data
    fitted_data <- data %>%
      dplyr::filter(!is.na(Leaf), !is.na(.data[[y_var]])) %>%
      dplyr::group_by(pot, genotype, treatment) %>%
      dplyr::filter(n() >= 4) %>%
      tidyr::nest() %>%
      
      dplyr::mutate(
        fit = purrr::map(data, ~ tryCatch(
          mgcv::gam(
            as.formula(paste0(y_var, " ~ s(Leaf, bs = 'cs', k = ", k, ")")),
            data = .x,
            method = "REML"
          ),
          error = function(e) NULL
        )),
        
        pred_data = purrr::map(data, ~
                                 grid %>%
                                 dplyr::filter(
                                   Leaf >= min(.x$Leaf),
                                   Leaf <= max(.x$Leaf)
                                 )
        ),
        
        y_pred = purrr::map2(fit, pred_data, ~ {
          if (is.null(.x)) return(NULL)
          pmax(0, as.numeric(predict(.x, newdata = .y)))
        }),
        
        fit_success = !purrr::map_lgl(fit, is.null)
      )
    
    # Print convergence summary for this variable
    cat("\n=== Convergence summary for", y_var, "===\n")
    print(table(fitted_data$fit_success))
    
    convergence_by_group <- fitted_data %>%
      dplyr::group_by(genotype, treatment) %>%
      dplyr::summarise(
        n_total = n(),
        n_failed = sum(!fit_success),
        n_success = sum(fit_success),
        .groups = "drop"
      )
    print(convergence_by_group)
    
    # Continue with filtering and returning results
    fitted_data %>%
      dplyr::filter(fit_success) %>%
      dplyr::select(pot, genotype, treatment, pred_data, y_pred) %>%
      tidyr::unnest(c(pred_data, y_pred)) %>%
      dplyr::mutate(variable = y_var)
  }
  
  purrr::map_dfr(
    y_vars,
    ~ create_single_pspline(data, .x, grid, k)
  )
}

# Usage
leaf_grid <- tibble(
  Leaf = seq(1, 22, by = 0.1))

#df
df_pspline_all <- create_multiple_pspline(
  data = my_dat0,
  y_vars = c("llength", "sheathl" ,"idist"),
  grid = leaf_grid,
  k = 11) %>%
dplyr::group_by(pot, variable) %>%
dplyr::mutate(norm_length=(y_pred - min(y_pred,na.rm = TRUE))/(max(y_pred,na.rm = TRUE) - min(y_pred,na.rm = TRUE)))

#adding relative to ear
df_pspline_all_ear<-df_pspline_all%>%
                    left_join(my_dat0%>%
                                dplyr::select(pot, earleaf)%>%
                                unique, by="pot")%>%
                    dplyr::group_by(pot)%>%
                    dplyr::mutate(rel_e=Leaf-earleaf)%>%
                    ungroup()%>%
  mutate(variable = recode(variable, !!!recode_map),
         variable = factor(variable, levels = panel_order))



# ==============================================================================
# 4. Combine with canopy dimensions
# ==============================================================================


recode_map <- c(llength = "Blade", sheathl = "Sheath", idist = "Internode")

panel_order <- c("Blade", "Sheath", "Internode")


mergePred <- int_phyllo %>%
  bind_rows(int_colla) %>%
  left_join(thermal_ti, by=c("pot", "genotype")) %>%
  dplyr::group_by(pot, parameter) %>%
  dplyr::mutate(intersection_point=case_when(intersection_point < 0~0,
                                             TRUE ~ intersection_point),
                intersection_point=case_when(parameter=="colla" & leaf<3| parameter=="phyllo" & leaf<3 ~NA,
                                             TRUE ~ intersection_point)) %>%
  ungroup()%>%
  dplyr::group_by(pot, leaf)%>%
  dplyr::mutate(flexposure=intersection_point-thrmal_intrsect) %>%
  ungroup() %>%
  #dplyr::rename(leaf)%>%
  left_join(df_pspline_all%>%
            dplyr::rename(leaf=Leaf)%>%
            dplyr::select(pot,leaf, y_pred,norm_length, variable), by=c("leaf", "pot"))%>%
  dplyr::filter(variable!="NA") %>%
  dplyr::mutate(
      type = dplyr::case_when(
      variable == "llength" ~ "phyllo",
      variable %in% c("sheathl", "idist") ~ "colla",
      TRUE ~ NA_character_), 
        dplyr::across(
          .cols = c(y_pred, norm_length),       # columns to set
          .fns = ~ ifelse(type != parameter, NA_real_, .x)))%>%
  dplyr::filter(!is.na(y_pred))%>%
  mutate(variable = recode(variable, !!!recode_map),
         variable = factor(variable, levels = panel_order))



# ==============================================================================
# 5.Check plot aspects
# ==============================================================================
distance <- mergePred %>%
  dplyr::filter(type == "phyllo" & variable == "llength") %>%
  dplyr::filter(norm_length == 1) %>%
  dplyr::summarize(
    exp_flex  = mean(flexposure, na.rm = TRUE),
    exp_range = IQR(flexposure, na.rm = TRUE),
    Q1        = quantile(flexposure, 0.25, na.rm = TRUE),
    Q3        = quantile(flexposure, 0.75, na.rm = TRUE),
    lower     = Q1 - 1.5 * IQR(flexposure, na.rm = TRUE),
    upper     = Q3 + 1.5 * IQR(flexposure, na.rm = TRUE)
  )
# ==============================================================================
# 6. Plots
# ==============================================================================

#separately by phytomer component
p5 <- ggplot(mergePred, aes(x = flexposure, y = y_pred, group = pot)) +
  geom_vline(xintercept = 0, linewidth = 1.2) +
  
  #individual plants
  geom_line(alpha = 0.3, colour = "grey60", linewidth = 0.5) +
  
  # Population smooths
  geom_smooth(
    data = \(x) subset(x, variable == "Blade"),   # Blade
    aes(group = 1),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    colour = "#E41A1C",
    linewidth = 1.2
  ) +
  
  geom_smooth(
    data = \(x) subset(x, variable == "Sheath"),   # Sheath
    aes(group = 1),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    colour = "#377EB8",
    linewidth = 1.2
  ) +
  
  geom_smooth(
    data = \(x) subset(x, variable == "Internode"),     # Internode
    aes(group = 1),
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = FALSE,
    colour = "#4DAF4A",
    linewidth = 1.2
  ) +
  
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  
  scale_x_continuous(
    limits = c(-25, 125),
    breaks = seq(-25, 125, by = 25),
    name = expression("Time from tassel initiation (d"[20*degree*C]*")")
  ) +
  
  
  # Custom y-axis per facet
  ggh4x::facetted_pos_scales(
    y = list(
      llength = scale_y_continuous(
        limits = c(0, 142),
        breaks = function(lim) seq(0, 140, 20)
      ),
      sheathl = scale_y_continuous(
        limits = c(0, 35),
        breaks = function(lim) seq(0, 35, 5)
      ),
      idist = scale_y_continuous(
        limits = c(0, 30),
        breaks = function(lim) seq(0, 30, 5)
      )
    )
  ) +
  
  labs(title = "A", y = "Length (cm)") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.minor = element_blank()
  )


#print it
p5

#save it
timestamp <- format(Sys.time(), "%d%m%Y")
ggsave(plot=p5, paste0(fig_dir, "Fig3A_", timestamp, ".pdf"), width = 7, height = 5, dpi = 300,  bg = "transparent")
saveRDS(p5, file = paste0(fig_dir, "Fig3A_", timestamp, ".rds"))


#Normalized plot combining all phytomers


p6 <- ggplot(mergePred %>% filter(norm_length >= 0),  # Cut at 0
       aes(x = flexposure, y = norm_length, 
           group = variable, 
           color = variable)) +  # Add color aesthetic
  
  geom_vline(xintercept = 0, linewidth = 1.2, colour = "black") +

  
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, bs = "cs"),
    se = TRUE,
    linewidth = 1.6
  ) +
  
  scale_colour_manual(
    values = c(
      "Blade"      = "#E41A1C",
      "Sheath"     = "#377EB8",
      "Internode"  = "#4DAF4A"
    ),
    labels = c(
      "Blade"       = "Blade length (normalized)",
      "Sheath"      = "Sheath length  (normalized)",
      "Internode"  = "Internode length (normalized)"
    ),
    name = "Phytomer component"  # Legend title
  ) +
  
  scale_x_continuous(
    limits = c(-25, 125),
    breaks = seq(-25, 125, by = 25),
    name = expression("Time from tassel initiation (d"[20*degree*C]*")")
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  
  # Add facet with your desired title
  facet_wrap(~ "Phytomer Component Profiles", strip.position = "top") +
  
  labs(
    title = "B",
    x = expression("Time from tassel initiation (d"[20*degree*C]*")"),
    y = "Normalized length"
  ) +
  
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    legend.title =  element_text(size=11, face="bold"),
    legend.text = element_text(size=12),
    plot.title = element_text(hjust = 0, size=15, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 90, vjust=0.5, hjust=1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size=14, margin = margin(t = 10)),
    axis.title.y = element_text(size=14, margin = margin(r = 10)),
    # panel.grid.minor = element_line(linewidth = 0.3, colour = "grey85")
    panel.grid.major.x = element_line(linewidth = 0.3, colour = "grey85"),  # Remove vertical major lines
    panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),  # Remove horizontal major lines
    panel.grid.minor.x = element_blank(),  # Remove vertical minor lines
    panel.grid.minor.y = element_blank()  # Remove horizontal minor lines
  )

#print it
p6

#save it 
timestamp <- format(Sys.time(), "%d%m%Y")
ggsave(plot=p6, paste0(fig_dir, "Fig3B_", timestamp, ".pdf"), width = 7/1.48, height = 5/1.48, dpi = 300,  bg = "transparent")
saveRDS(p6, file = paste0(fig_dir, "Fig3B_", timestamp, ".rds"))


