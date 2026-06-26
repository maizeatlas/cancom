rm(list = ls()) #clean environment
gc() #garbage collection

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# Leaf-tip formalism benchmarking: empirical vs. predicted T.I. timing
#
#   Pipeline:
#     1. SAM dissection data (empirical T.I.) and final leaf number BLUEs are
#        merged and six published leaf-counting formalisms are applied
#        (Johnson, Lejeune, Padilla, Drouault, Tollenaar, Arnold hi/lo).
#     2. Each formalism is benchmarked against empirical T.I. via linear
#        regression, yielding R², adjusted R², intercept, slope, RMSE, MAE,
#        and bias per model.
#     3. The best-fitting formalism (Johnson) is plotted as observed vs.
#        predicted, coloured by treatment (control / night-break), and
#        assembled with a companion panel (Fig1A) into a final two-panel figure.
#
#   Outputs:
#     A. Table of model performance statistics across all formalisms.
#     B. Fig2B — observed vs. predicted scatter for the Johnson formalism.
#     C. Fig2 — two-panel composite (Fig2A | Fig2B) saved as PDF.
# ==============================================================================

#packages
library(tidyverse) #data manipulation and ggplot
library(lme4) #mixed model
library(knitr) #table
library(ggrepel) #plotting
library(patchwork) #combining

# ---- Paths ----
dat_dir = paste0(getwd(),"/data/")
fig_dir = paste0(getwd(),"/figures/")

# ==============================================================================
# 1. Read and summarize data 
# ==============================================================================

#leaf tip at which T.I. occurred from dissections
samint_alt <- read.csv(paste0(dat_dir,"02_dat_platform_SAM.csv")) %>%  
  dplyr::mutate(treatment=as.factor(fct_recode(treatment, "control"="con", "nightbreak"="nb"))) %>%
  rename(empirical_ti=fit)

#total leaf number BLUEs from harvest
tot_leaf <- read.csv(paste0(dat_dir,"02_dat_platform_finalleafno.csv")) %>%
  dplyr::mutate(across(c(genotype,treatment),as.factor),
                treatment=fct_recode(treatment, "nightbreak"="NB", "control"="Control"))

#merge and add formalisms
mrgelftran<-samint_alt%>%
  left_join(tot_leaf, by=c("genotype", "treatment"))%>%
  dplyr::mutate(Johnson=(emmean-6.5)/1.51, #This study
                Lejeune=(emmean-1.95)/1.84,  #Lejeune & bernier
                Padilla=(-16.09+(2.36*8.04))+((emmean-8.04)*0.63),#Padilla & Otegui
                Drouault=(emmean-5.5)/1.73, #Drouault
                Tollenaar=0.5*emmean, #Tollenaar
                Arnold_hi=0.54*emmean, #Arnold (hi)
                Arnold_lo=0.42*emmean)%>% #Arnold  (lo)
  droplevels()%>%
  dplyr::filter(genotype !="CML258")


# ==============================================================================
# 2. Compare formalisms
# ==============================================================================

stats_by_form_lm <- mrgelftran %>%
  dplyr::select(genotype, empirical_ti, treatment, emmean, Johnson:Arnold_lo) %>%
  pivot_longer(c(Johnson, Lejeune:Arnold_lo), values_to = "predict_ti", names_to = "mod") %>%
  dplyr::group_by(mod) %>%
  dplyr::reframe({
    model <- lm(predict_ti ~ empirical_ti, na.action = na.exclude)
    model_summary <- summary(model)
    
    tibble(
      R2 = model_summary$r.squared,
      R2_adj = model_summary$adj.r.squared,
      Intercept = coef(model)[1],
      Slope = coef(model)[2],
      p_value = model_summary$coefficients[2, 4],
      RMSE = sqrt(mean((predict_ti - empirical_ti)^2, na.rm = TRUE)),
      MAE = mean(abs(predict_ti - empirical_ti), na.rm = TRUE),
      Bias = mean(predict_ti - empirical_ti, na.rm = TRUE)
    )
  })

#print table
stats_by_form_lm %>%
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  kable(
    col.names = c("Model", "R²", "R2_adj", "Intercept", "Slope", "p_value", "RMSE", "MAE", "Bias"),
    caption = "Model Performance Statistics by Formulation",
    align = c("l", "c", "c", "c", "c", "c", "c", "c", "c"),
    format = "pipe"  # or "html", "latex"
  )

# ==============================================================================
# 3. Plot
# ==============================================================================

# Extract values for Johnson estimates
johnson_stats <- stats_by_form_lm %>% filter(mod == "Johnson")

# Calculate the common axis range based on the maximum value
axis_max <- max(c(mrgelftran$Johnson, mrgelftran$empirical_ti), na.rm = TRUE)
axis_min <- min(c(mrgelftran$Johnson, mrgelftran$empirical_ti), na.rm = TRUE)

p2 <- ggplot(mrgelftran,
       aes(x = empirical_ti, y = Johnson, color = treatment)) +
  scale_x_continuous(
    limits = c(axis_min, axis_max),
    minor_breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1)
  ) +
  scale_y_continuous(
    limits = c(axis_min, axis_max),
    minor_breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1)
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    colour = "black",
    linewidth = 1,
    alpha = 0.5,
    linetype = "dashed"
  ) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_discrete(
    name = "Treatment",
    labels = c(
      "control" = "Control",
      "nightbreak" = "Night-Break"
    )
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    colour = "black",
    linewidth = 1.2
  ) +
  labs(
    title = "B",
    x = "Observed (visible leaf at T.I.)",
    y = "Predicted (visible leaf at T.I.)",
    color = "Treatment"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.28, vjust = 1.8,
    size = 5,
    parse = TRUE,
    label = sprintf('italic(R)^2~"="~"%.2f"', johnson_stats$R2_adj)
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.21, vjust = 4.3,
    size = 5,
    parse = TRUE,
    label = sprintf('RMSE~"="~"%.2f"', johnson_stats$RMSE)
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    legend.title =  element_text(size=11, face="bold"),
    legend.text = element_text(size=12),
    legend.box = "vertical",
    legend.box.just = "top",
    legend.margin = margin(t = -155),
    plot.title = element_text(hjust = 0, size = 15, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    panel.grid.minor = element_line(linewidth = 0.3, colour = "grey85")
  )
p2

timestamp <- format(Sys.time(), "%d%m%Y")
ggsave(plot=p2, paste0(fig_dir, "Fig2B_", timestamp, ".pdf"), width = 7, height = 5, dpi = 300)
saveRDS(p2, file = paste0(fig_dir, "Fig2B_", timestamp, ".rds"))


# install.packages("Polychrome")
# install.packages("ggrepel")
### to pull later for assembling plots

# ==============================================================================
# 4. Combine
# ==============================================================================

p1 <- readRDS(paste0(fig_dir, "Fig2A_01062026", ".rds"))
p2 <- readRDS(paste0(fig_dir, "Fig2B_18062026", ".rds"))
final_plot <- (p1 | plot_spacer() | p2) +
  plot_layout(widths = c(1, 0.06, 1), guides = "keep")
final_plot <- p1 + 
  theme(plot.margin = margin(r = 15)) | 
  p2
final_plot
timestamp <- format(Sys.time(), "%d%m%Y")
ggsave(plot=final_plot, paste0(fig_dir, "Fig2_", timestamp, ".pdf"), width = 11, height = 5, dpi = 300)



