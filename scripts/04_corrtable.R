rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# Phytomer-curve landmark extraction: blade/sheath/internode growth curves and
# their relationship to the vegetative-to-reproductive transition (t_i)
#
#   Pipeline:
#     1. Platform (pot) BLUEs are read in, t_i (vegetative-to-reproductive
#        transition index, from tot_leaf) and rel_t are derived, genotype
#        CML258 is dropped, and per-pot observation counts for blade length
#        (llength), sheath length (sheathl), and internode distance (idist)
#        are tabulated (mean/min/max per trait).
#     2. Per-pot growth curves are fit against Leaf rank for each of the three
#        traits using statgenHTP::detectSingleOut (locfit-based smoothing,
#        nnLocfit = 0.8), yielding predicted curves (yPred); the three trait
#        results are merged into one long-format data frame
#        (mrge_splne_cnpy).
#     3. Sheath plateau detection: for each pot, the sheath-length curve's
#        post-maximum decline is isolated, and changepoint::cpt.meanvar
#        (BinSeg, manual penalty) locates where elongation decline stops,
#        giving a per-pot "shthplat" leaf rank (falls back to the
#        smallest-decline point if no changepoint is found).
#     4. Per-pot/trait summary stats (max/min value and their Leaf-rank
#        position) are extracted from the fitted curves and joined with the
#        sheath plateau, t_i, earleaf, frst_idist_rnk, and genotype.
#     5. Correlation analysis: Pearson correlations (with significance
#        stars) between curve-derived timing features (max_x position,
#        frst_idist_rnk, shthplat) and t_i / earleaf are computed per trait
#        (blade/sheath/internode), with biologically irrelevant
#        feature × trait combinations filtered out (e.g., frst_idist_rnk vs.
#        sheath/blade, shthplat vs. internode/blade).
#
#   Outputs:
#     A. Console table — mean/min/max number of non-NA observations per pot
#        for llength, sheathl, idist.
#     B. Console table (cor_table_tidy) — correlation coefficients,
#        p-values, and significance codes linking curve-timing features
#        (max_x, frst_idist_rnk, shthplat) to t_i and earleaf, by trait.
#     (No files are written to disk in this script — all results are
#      printed for inspection; a commented-out block at the end sketches an
#      alternative raw-data correlation approach but is not executed.)
# ==============================================================================

library(tidyverse) #data manipulation
library(purrr)
library(mgcv)
library(statgenHTP)
library(segmented)
library(changepoint)

# ---- Paths ----
dat_dir = paste0(getwd(),"/data/")
fig_dir = paste0(getwd(),"/figures/")


# ==============================================================================
# 1. Read and summarize data 
# ==============================================================================

# Modeled unadjusted BLUEs from different experiments
my_dat0 <- read.csv(paste0(dat_dir,"04_dat_platform_canopy.csv"))%>%
           dplyr::mutate(t_i=(tot_leaf-6.5)/1.51)%>%
          dplyr::group_by(pot)%>%
          dplyr::mutate(rel_t=Leaf-t_i)%>%
           dplyr::mutate(pot=as.factor(pot))%>%
  dplyr::filter(genotype!="CML258")


#summary
cat("Mean number of observations per pot:\n\n")

# Helper to count non-NA observations per pot
count_per_pot <- function(var) {
  aggregate(
    !is.na(var) ~ pot,
    data = my_dat0,
    FUN = sum
  )
}

# Counts per pot for each variable
llength_pot  <- count_per_pot(my_dat0$llength)
sheathl_pot  <- count_per_pot(my_dat0$sheathl)
idist_pot    <- count_per_pot(my_dat0$idist)

# Mean across pots
mean_obs <- data.frame(
  Variable = c("llength", "sheathl", "idist"),
  Mean_observations_per_pot = c(
    mean(llength_pot[, 2]),
    mean(sheathl_pot[, 2]),
    mean(idist_pot[, 2])
  )
)

print(mean_obs)
cat("\n")

obs_range <- data.frame(
  Variable = c("llength", "sheathl", "idist"),
  Mean = c(mean(llength_pot[,2]), mean(sheathl_pot[,2]), mean(idist_pot[,2])),
  Min  = c(min(llength_pot[,2]), min(sheathl_pot[,2]), min(idist_pot[,2])),
  Max  = c(max(llength_pot[,2]), max(sheathl_pot[,2]), max(idist_pot[,2]))
)

print(obs_range)

# ==============================================================================
# 2.Phytomer curves
# ==============================================================================

#create timepoint, so that it can be passed on to fit
TPcanopy <- createTimePoints(dat = my_dat0,
                            experimentName = "Phenoarch",
                            genotype = "genotype",
                            timePoint = "Leaf",
                            #repId = "rep",
                            plotId = "pot")
                            #rowNum = "Line",
                            #colNum = "NumParRang")

#blade
resublade <- detectSingleOut(TP = TPcanopy,
                             trait = "llength",
                             #plotIds = c("1473"),
                             confIntSize = 5,
                             nnLocfit = 0.8)%>%
  dplyr::mutate(timestamp = ymd_hms(timePoint),      
                Leaf = second(timePoint), 
                variable="blade") 

#sheath
resusheath <- detectSingleOut(TP = TPcanopy,
                              trait = "sheathl",
                              #plotIds = c("1473"),
                              confIntSize = 5,
                              nnLocfit = 0.8)%>%
  dplyr::mutate(timestamp = ymd_hms(timePoint),      
                      Leaf = second(timePoint), 
                variable="sheath")




#internode
resuinter <- detectSingleOut(TP = TPcanopy,
                              trait = "idist",
                              #plotIds = c("1473"),
                              confIntSize = 5,
                              nnLocfit = 0.8)%>%
  dplyr::mutate(timestamp = ymd_hms(timePoint),      
                Leaf = second(timePoint),
                variable="inter") 

#merge
mrge_splne_cnpy<-  resusheath %>%
                   bind_rows(resublade)%>%
                   bind_rows(resuinter)%>%
                    dplyr::select(
                      -contains("time"),        # remove columns that contain "time"
                      -sheathl, -llength, -idist, -outlier)%>%# remove these specific columns
                    dplyr::rename(pot=plotId)


# ==============================================================================
# 2. Extracting sheath breakpoint (plateau)
# ==============================================================================

#plot to visualize sheaths after max
# resusheath %>%
#   filter(plotId %in% unique(resusheath$plotId)[1:9]) %>%
#   group_by(plotId) %>%
#   mutate(max_leaf = Leaf[which.max(yPred)]) %>%
#   ggplot(aes(x = Leaf, y = yPred)) +
#   geom_line() +
#   geom_point(size = 1) +
#   geom_vline(aes(xintercept = max_leaf), color = "red", linetype = "dashed") +
#   facet_wrap(~plotId, scales = "free") +
#   theme_minimal() +
#   labs(title = "yPred vs Leaf (red line = maximum)")


# xyplot(yPred~Leaf|plotId, data=resusheath)+
#   layer(panel.abline(v=c(5,10,15)))#, subset=plotId=="1669")


fit_cpt_shth <- function(df) {
  tryCatch({
    # Sort by Leaf
    df <- df %>% arrange(Leaf)
    
    # Find maximum
    max_idx <- which.max(df$yPred)
    
    # Need at least 5 points after maximum
    if (max_idx >= nrow(df) - 4) {
      return(NULL)
    }
    
    # Get post-maximum data
    post_max <- df[max_idx:nrow(df), ]
    
    # Calculate the rate of change (first derivative approximation)
    diff_y <- -diff(post_max$yPred)  # negative because we expect decline
    
    # Use changepoint on the rate of change to find where decline stops
    if (length(diff_y) < 3) return(NULL)
    
    cpt_mod <- cpt.meanvar(diff_y, 
                           method = "BinSeg",
                           Q = 1,  # Look for 1 changepoint
                           penalty = "Manual",
                           pen.value = 0.05)
    
    # Store data for later extraction
    attr(cpt_mod, "post_max_df") <- post_max
    attr(cpt_mod, "max_idx") <- max_idx
    
    return(cpt_mod)
    
  }, error = function(e) {
    message("Error in plotId: ", unique(df$plotId), " - ", e$message)
    return(NULL)
  })
}

# Fit models
segmodelshth <- resusheath %>% 
  dplyr::group_by(plotId) %>%
  arrange(Leaf) %>%
  group_modify(~ tibble(model = list(fit_cpt_shth(.))))

# Extract breakpoint
extract_brkpointsth <- function(model) {
  if (is.null(model)) {
    return(NA)
  }
  
  tryCatch({
    cpts_idx <- cpts(model)
    post_max_df <- attr(model, "post_max_df")
    
    if (length(cpts_idx) == 0) {
      # If no changepoint detected, use the point where decline is smallest
      diff_y <- -diff(post_max_df$yPred)
      plateau_idx <- which.min(abs(diff_y)) + 1
    } else {
      # Add 1 because diff reduces length by 1
      plateau_idx <- cpts_idx[1] + 1
    }
    
    # Make sure we don't go out of bounds
    plateau_idx <- min(plateau_idx, nrow(post_max_df))
    
    return(post_max_df$Leaf[plateau_idx])
    
  }, error = function(e) {
    return(NA)
  })
}

# Create results dataframe
breakvsth_df <- segmodelshth %>%
  dplyr::mutate(
    shthplat = map_dbl(model, extract_brkpointsth),
    pot = as.factor(plotId)
  ) %>%
  ungroup()%>%
  dplyr::select(pot, shthplat)

# Check results

#print(breakvsth_df)
#summary(breakvsth_df$shthplat)



# ==============================================================================
# 3. Max, min leaf blade, sheath, and internode
# ==============================================================================

# Extract summary statistics 
summary_curve <- mrge_splne_cnpy %>%
  dplyr::group_by(pot, variable) %>%
  dplyr::summarise(
    max_value = max(yPred, na.rm = TRUE),
    min_value = min(yPred, na.rm = TRUE),
    max_x = Leaf[which.max(yPred)],   # x corresponding to max y
    min_x = Leaf[which.min(yPred)],   # x corresponding to min y
    .groups = "drop")%>%  
   left_join(breakvsth_df, by="pot")%>%
   left_join(my_dat0 %>%
              dplyr::select(t_i,earleaf,frst_idist_rnk, genotype, pot)%>%
              distinct(pot, .keep_all = TRUE), by = "pot")

# ==============================================================================
# 4. Correlation table+significance
# ==============================================================================

cor_table_tidy <- summary_curve %>%
  pivot_longer(
    cols = c("max_value":"shthplat", "frst_idist_rnk"),
    names_to = "feature",
    values_to = "value")  %>%
  dplyr::filter(feature %in% c("max_x","frst_idist_rnk", "shthplat"))%>%
  dplyr::group_by(variable, feature) %>%
  dplyr::summarise(
    #n = n(),
    cor_ti = cor(value, t_i, use = "complete.obs"),
    p_ti = cor.test(value, t_i)$p.value,
    cor_ear = cor(value, earleaf, use = "complete.obs"),
    p_ear   = cor.test(value, earleaf)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    sig_ti = case_when(
      p_ti < 0.001 ~ "***",
      p_ti < 0.01  ~ "**",
      p_ti < 0.05  ~ "*",
      TRUE ~ "ns"
    ),
    sig_ear = case_when(
      p_ear < 0.001 ~ "***",
      p_ear < 0.01  ~ "**",
      p_ear < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )%>%
     dplyr::filter(
    !(feature == "frst_idist_rnk" &
        variable %in% c("sheath", "blade")))%>%
      dplyr::filter(
          !(feature == "shthplat" &
              variable %in% c("inter", "blade")))%>%
  dplyr::select(variable, feature, cor_ti, p_ti, sig_ti, cor_ear, p_ear, sig_ear )
      
#print
print(cor_table_tidy)


#raw data correlations

# summary_curve <- my_dat0 %>%
#   pivot_longer(
#     cols = c("idist":"llength"),
#     names_to = "feature",
#     values_to = "value") %>%
#   dplyr::group_by(pot,genotype, feature) %>%
#   dplyr::summarise(
#     max_value = max(value, na.rm = TRUE),
#     min_value = min(value, na.rm = TRUE),
#     max_x = Leaf[which.max(value)],   # x corresponding to max y
#     min_x = Leaf[which.min(value)],   # x corresponding to min y
#     .groups = "drop")%>%
#   left_join(my_dat0 %>%
#               dplyr::select(t_i,earleaf, pot)%>%
#               distinct(pot, .keep_all = TRUE), by = "pot")%>%
#   dplyr::group_by(feature) %>%
#   dplyr::summarise(
#     n = n(),
#     cor_ti = cor(max_x, t_i, use = "complete.obs"),
#     p_ti = cor.test(max_x, t_i)$p.value,
#     cor_ear = cor(max_x, earleaf, use = "complete.obs"),
#     p_ear   = cor.test(max_x, earleaf)$p.value,
#     .groups = "drop"
#   )

