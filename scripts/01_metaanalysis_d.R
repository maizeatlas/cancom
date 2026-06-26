rm(list = ls())
gc()

# Set the working directory
setwd("/rootpath/to/your/project")

# ==============================================================================
# Plastochron–phyllochron scaling: meta-analytic estimation of the plasto/phyllo
# slope across published maize studies
#
#   Pipeline:
#     1. Per-study linear models regress leaf primordia number (plasto) on leaf
#        tip number (phyllo), with genotype included as a covariate in
#        multi-genotype studies; slopes and standard errors are extracted for
#        each study.
#     2. A fixed-effects meta-analysis pools slopes via inverse-variance
#        weighting, yielding a pooled slope and 95% CI across all studies.
#     3. Two additional models benchmark the pooled estimate: a simple OLS
#        controlling for study (origin) and a nested OLS with genotype nested
#        within study; AIC and R² are compared across all three approaches.
#     4. Genotype-level slope heterogeneity is further probed within the Padilla and Otegui
#        study using a phyllo × genotype interaction model and emmeans-derived
#        per-genotype trends (both unweighted and inverse-variance weighted).
#
#   Outputs:
#     A. Console summaries of per-study slopes, pooled meta-analytic estimate,
#        nested/simple model fits, and genotype-level trend estimates.
#     B. Fig2A — plasto vs. phyllo scatter coloured by study (point opacity
#        scaled by inverse-variance weight) with the pooled meta-analytic
#        regression line and annotated slope ± SE.
# ==============================================================================

library(tidyverse) #data manipulation and ggplot
library(data.table)

# ---- Paths ----
dat_dir = paste0(getwd(),"/data/")
fig_dir = paste0(getwd(),"/figures/")

# ==============================================================================
# 1. Read and summarize data 
# ==============================================================================

# Modeled unadjusted BLUEs from different experiments
my_dat0 = fread(file = paste0(dat_dir, "01_dat_metanalaysis.csv"))

my_dat1 <- my_dat0 %>%
  filter(!origin %in% c("Thiagarajah", "Warrington", "Zur"))

# Meta-analysis of plasto ~ phyllo relationship across studies
# Assumes my_dat1 has columns: origin (study), genotype, phyllo (x), plasto (y)

# Check data structure
cat("Data Summary:\n")
cat("Total observations:", nrow(my_dat1), "\n")
cat("Number of studies:", length(unique(my_dat1$origin)), "\n")
cat("Studies:", paste(unique(my_dat1$origin), collapse = ", "), "\n\n")

# Summary by study
study_summary <- aggregate(cbind(n = plasto) ~ origin, data = my_dat1, FUN = length)
geno_summary <- aggregate(cbind(n_genotypes = genotype) ~ origin, data = my_dat1, 
                          FUN = function(x) length(unique(x)))
study_info <- merge(study_summary, geno_summary)
colnames(study_info) <- c("Study", "N_observations", "N_genotypes")
print(study_info)
cat("\n")

# ==============================================================================
# 2. Meta-analysis type I (inverse-variance weighting)
# ==============================================================================

# Fit individual linear models for each study
# For studies with multiple genotypes, account for genotype effects
studies <- unique(my_dat1$origin)
models_list <- list()
slopes <- numeric(length(studies))
se_slopes <- numeric(length(studies))

for(i in seq_along(studies)) {
  study_data <- subset(my_dat1, origin == studies[i])
  n_genotypes <- length(unique(study_data$genotype))
  
  # If multiple genotypes, include genotype as covariate
  if(n_genotypes > 1) {
    models_list[[i]] <- lm(plasto ~ phyllo + genotype, data = study_data)
  } else {
    models_list[[i]] <- lm(plasto ~ phyllo, data = study_data)
  }
  
  slopes[i] <- coef(models_list[[i]])[2]  # plasto coefficient
  se_slopes[i] <- summary(models_list[[i]])$coef[2, 2]
}

names(models_list) <- studies

cat("Individual Study Results (accounting for genotype when applicable):\n")
for(i in seq_along(studies)) {
  study_data <- subset(my_dat1, origin == studies[i])
  n_genotypes <- length(unique(study_data$genotype))
  geno_note <- ifelse(n_genotypes > 1, 
                      paste0(" (", n_genotypes, " genotypes)"), 
                      "")
  cat(studies[i], geno_note, "slope:", round(slopes[i], 3), 
      "±", round(se_slopes[i], 3), "\n")
}

# Meta-analysis: Fixed-effects model (inverse-variance weighted)
weights <- 1 / se_slopes^2
pooled_slope <- sum(weights * slopes) / sum(weights)
pooled_se <- sqrt(1 / sum(weights))

cat("Meta-Analysis Results (Fixed Effects):\n")
cat("Pooled slope:", round(pooled_slope, 3), "\n")
cat("Standard error:", round(pooled_se, 3), "\n")
cat("95% CI: [", round(pooled_slope - 1.96*pooled_se, 3), ",", 
    round(pooled_slope + 1.96*pooled_se, 3), "]\n\n")

# ==============================================================================
# 2. Meta-analysis type II (nested and simple linear models)
# ==============================================================================

# Combined models with genotype nested within study
# Create a unique genotype identifier nested within study
my_dat1$study_geno <- paste(my_dat1$origin, my_dat1$genotype, sep = "_")

nested_model <- lm(plasto ~ phyllo + study_geno, data = my_dat1)
cat("Nested Model (genotype nested within study):\n")
print(summary(nested_model))

# Also fit model with just study (for comparison)
simple_model <- lm(plasto ~ phyllo + origin, data = my_dat1)
cat("\n\nSimple Model (controlling for study only):\n")
print(summary(simple_model))

# Compare models
cat("\n\nModel Comparison:\n")
cat("Nested model AIC:", AIC(nested_model), "\n")
cat("Simple model AIC:", AIC(simple_model), "\n")
cat("Nested model R-squared:", summary(nested_model)$r.squared, "\n")
cat("Simple model R-squared:", summary(simple_model)$r.squared, "\n\n")

# Test for genotype effects using Padilla study
my_dat1_Padilla <- my_dat1 %>%
  filter(origin == "Padilla")
my_dat1_Padilla$genotype <- factor(my_dat1_Padilla$genotype)

geno_model <- lm(plasto ~ phyllo + genotype + phyllo:genotype, data = my_dat1_Padilla)
summary(geno_model)

library(emmeans)
phyllo_trends <- emtrends(geno_model, ~ genotype, var = "phyllo")
phyllo_trends_df <- as.data.frame(phyllo_trends)

mean_slope <- mean(phyllo_trends_df$phyllo.trend)
mean_slope
se_mean <- sd(phyllo_trends_df$phyllo.trend) / sqrt(nrow(phyllo_trends_df))
se_mean

# weights = 1 / SE^2
w <- 1 / (phyllo_trends_df$SE^2)
weighted_mean <- sum(w * phyllo_trends_df$phyllo.trend) / sum(w)
weighted_se <- sqrt(1 / sum(w))

weighted_mean
weighted_se

# ==============================================================================
# 3. Plot results
# ==============================================================================

# Create data frame for regression lines
pred_range <- seq(min(my_dat1$phyllo), max(my_dat1$phyllo), length.out = 100)

# 1. Meta-analysis line
intercept_meta <- mean(my_dat1$plasto) - pooled_slope * mean(my_dat1$phyllo)
meta_line <- data.frame(
  phyllo = pred_range,
  plasto = intercept_meta + pooled_slope * pred_range,
  model = "Meta-analysis (pooled)"
)

# 2. Simple model line
simple_pred <- predict(simple_model, 
                       newdata = data.frame(phyllo = pred_range, 
                                            origin = studies[1]))
simple_line <- data.frame(
  phyllo = pred_range,
  plasto = simple_pred,
  model = "Simple model (study)"
)

# 3. Nested model line
nested_pred <- predict(nested_model, 
                       newdata = data.frame(phyllo = pred_range, 
                                            origin = studies[1],
                                            genotype = unique(my_dat1$genotype)[1],
                                            study_geno = paste(studies[1], 
                                                               unique(my_dat1$genotype)[1], 
                                                               sep = "_")))
nested_line <- data.frame(
  phyllo = pred_range,
  plasto = nested_pred,
  model = "Nested model (study+genotype)"
)

# Combine all lines
all_lines <- rbind(meta_line, simple_line, nested_line)
all_lines$model <- factor(all_lines$model, 
                          levels = c("Meta-analysis (pooled)", 
                                     "Simple model (study)", 
                                     "Nested model (study+genotype)"))

# Create the plot
p <- ggplot(my_dat1, aes(x = phyllo, y = plasto, color = origin)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_line(data = all_lines, 
            aes(x = phyllo, y = plasto, linetype = model, group = model),
            color = "black", linewidth = 1.2) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  labs(title = "Meta-Analysis: All Studies with Three Modeling Approaches",
       x = "Leaf tip number",
       y = "Leaf primordia number",
       color = "Study",
       linetype = "Model") +
  theme_bw() +
  theme(legend.position = "right",
        legend.box = "vertical",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=14, margin = margin(b = 10)),
        axis.title.y = element_text(size=14, margin = margin(r = 10)))

print(p)

# Create the manuscript plot
all_lines_d <- all_lines %>%
  filter(model == "Meta-analysis (pooled)")

# install.packages("showtext")
# library(showtext)
# font_add_google("Roboto Mono", "robomono")
# showtext_auto()

study_weights <- data.frame(
  origin = names(models_list),
  weight = weights / sum(weights)
)

study_weights$alpha <- scales::rescale(
  study_weights$weight,
  to = c(0.2, 1)  # adjust if you want more/less contrast
)

my_dat1_p <- my_dat1 %>%
  left_join(study_weights, by = "origin")

alpha_vals <- my_dat1_p %>%
  distinct(origin, alpha) %>%
  arrange(origin) %>%
  { setNames(.$alpha, .$origin) }

p1 <- ggplot(my_dat1_p, aes(phyllo, plasto, color = origin)) +
  geom_point(aes(alpha = origin), size = 2) +
  scale_color_discrete(name = "Study") +
  scale_alpha_manual(name = "Study", values = alpha_vals) +
  # guides(colour = guide_legend(override.aes = list(size = 4))) +
  geom_line(
    data = all_lines_d,
    aes(phyllo, plasto),
    color = "black",
    linewidth = 1.2
  ) +
  scale_x_continuous(
    minor_breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1)
  ) +
  scale_y_continuous(
    minor_breaks = function(limits) seq(floor(limits[1]), ceiling(limits[2]), by = 1)
  ) +
  labs(
    title = "A",
    x = "Leaf tip number",
    y = "Leaf primordia number",
    color = "Study"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    hjust = -0.15, vjust = 1.8,
    size = 5,
    parse = TRUE,
    label = sprintf(
      "hat(beta) == %.2f*'; '~SE(hat(beta)) == %.2f",
      pooled_slope, pooled_se
    )
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    legend.title =  element_text(size=11, face="bold"),
    legend.text = element_text(size=12),
    plot.title = element_text(hjust = 0, size=15, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.x = element_text(size=14, margin = margin(t = 10)),
    axis.title.y = element_text(size=14, margin = margin(r = 10)),
    panel.grid.minor = element_line(linewidth = 0.3, colour = "grey85")
  )
p1

timestamp <- format(Sys.time(), "%d%m%Y")
ggsave(plot=p1, paste0(fig_dir, "Fig2A_", timestamp, ".pdf"), width = 7, height = 5, dpi = 300)
saveRDS(p1, file = paste0(fig_dir, "Fig2A_", timestamp, ".rds"))

