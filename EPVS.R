library(readxl)
library(dplyr)
library(ggplot2)
library(boot)
library(patchwork)
library(lme4)
library(lmtest)
library(broom.mixed)
library(car)
library(MatchIt)
library(mgcv)
library(gamm4)
library(broom)

# Read in the dataset
scattering <- read_excel("~/Downloads/Datasets/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "scattering")
retardance <- read_excel("~/Downloads/Datasets/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "retardance")
orientation <- read_excel("~/Downloads/Datasets/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "orientation")


scattering$property <- "scattering"
retardance$property <- "retardance"
orientation$property <- "orientation"

main_data <- bind_rows(scattering, retardance, orientation)

# Clean variable names for plotting
main_data$Groups <- factor(main_data$Groups, levels = c("control", "experimental"))
main_data$Region <- factor(main_data$Region, levels = c("front", "occipital"))

# Compute block-level means

block_means <- main_data %>%
  group_by(Groups, Region, subID, property, distance) %>%
  summarise(mean_value = mean(OpticalProperty, na.rm = TRUE), .groups = "drop")


# EPVS are experimental adn non-healthy
# Vessel are controled
# Distance in microns

# fit linearity mixed model
# Continue EDA and models. 
# There are more observations in control than experiemental
# There are also more frontal thn occiptial
# Bootstrap amy be used
# Each distance should be model each coefficients
# Test for significance at each distance for each property
# Perform a likelihood ratio test
# Linear mixed model
# EDA

# Counts per group and region
main_data %>% 
  count(Groups, Region) %>% 
  arrange(desc(n)) %>% 
  arrange(desc(n))

# Counts per property and group
main_data %>% 
  count(property, Groups) %>% 
  arrange(desc(n))

# Some subjects may be missing some patients

ggplot(main_data, aes(x = OpticalProperty, fill = Groups)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  facet_wrap(~property + Region, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Optical Property by Group & Region",
       x = "Optical Property Value",
       y = "Count")


ggplot(main_data, aes(x = distance, y = OpticalProperty, color = Groups)) +
  geom_point(alpha = 0.4) +
  geom_smooth(se = FALSE) +
  facet_wrap(~property + Region, scales = "free_y") +
  theme_minimal() +
  labs(title = "Optical Property vs Distance",
       y = "Optical Property")

# Create summmary (mean +- SE) by group, region, and distance
summary_df <- summary_df %>%
  mutate(
    Region = case_when(
      is.na(Region) ~ "occipital",
      TRUE ~ as.character(Region)
    ),
    Region = factor(Region, levels = c("front", "occipital"))
  )
# Define plot functions 
plot_property <- function(df, region_filter, property_label) {
  df %>%
    filter(property == property_label, Region %in% region_filter) %>%
    ggplot(aes(x = distance, y = mean_val, color = Groups)) +
    geom_point(size = 2) +  # keep dots
    geom_errorbar(aes(ymin = mean_val - se, ymax = mean_val + se), width = 10) +
    # optional smooth trend line (remove if you want only dots + bars)
    # geom_smooth(method = "loess", se = FALSE, linetype = "dashed", alpha = 0.4) +
    scale_color_manual(values = c("black", "red")) +
    theme_minimal(base_size = 12) +
    labs(
      title = paste(paste(region_filter, collapse = " & "), property_label, "vs. Distance"),
      x = "Distance (µm)",
      y = property_label
    )
}



# Combined (all regions)
combined_summary <- summary_df %>%
  group_by(Groups, distance, property) %>%
  summarise(
    mean_val = mean(mean_val, na.rm = TRUE),
    se = sqrt(sum(se^2, na.rm = TRUE)) / 2,  # conservative SE estimate
    .groups = "drop"
  ) %>%
  mutate(Region = "combined")

summary_df_balanced <- bind_rows(summary_df, combined_summary)


p_combined_scatt <- plot_property(summary_df_balanced, "combined", "scattering")
p_combined_ret   <- plot_property(summary_df_balanced, "combined", "retardance")
p_combined_ori   <- plot_property(summary_df_balanced, "combined", "orientation")


p_front_scatt <- plot_property(summary_df, "front", "scattering")
p_front_ret   <- plot_property(summary_df, "front", "retardance")
p_front_ori   <- plot_property(summary_df, "front", "orientation")

p_occ_scatt <- plot_property(summary_df, "occipital", "scattering")
p_occ_ret   <- plot_property(summary_df, "occipital", "retardance")
p_occ_ori   <- plot_property(summary_df, "occipital", "orientation")

top_row <- p_combined_scatt + p_combined_ret + p_combined_ori
mid_row <- p_front_scatt + p_front_ret + p_front_ori
bottom_row <- p_occ_scatt + p_occ_ret + p_occ_ori

final_plot <- (top_row) / (mid_row) / (bottom_row) +
  plot_annotation(title = "Optical Properties vs Distance")

final_plot

#- ------ MIXED MODELS
set.seed(1234)

# Ensure factors / types
# --- Fix and standardize variable types / labels ---

main_data <- main_data %>%
  mutate(
    Groups = factor(Groups),             # control / experimental
    subID = factor(subID),
    distance = as.numeric(distance),
    OpticalProperty = as.numeric(OpticalProperty),
    property = factor(property, levels = c("scattering", "retardance", "orientation")),
    Region = tolower(as.character(Region))
  )

# Replace NA or blank with occipital, and recode synonyms properly
main_data <- main_data %>%
  mutate(
    Region = case_when(
      is.na(Region) ~ "occipital",
      Region %in% c("", "na") ~ "occipital",
      Region %in% c("front", "frontal") ~ "front",
      Region %in% c("occ", "occipital") ~ "occipital",
      TRUE ~ Region
    ),
    Region = factor(Region, levels = c("front", "occipital"))
  )

# Normalize Region factor to two levels: front and occipital
main_data <- main_data %>% mutate(Region = factor(ifelse(grepl("front", tolower(Region)), "front", "occipital"),
                                        levels = c("front","occipital")))

# Quick checks
cat("Unique properties:", paste(unique(main_data$property), collapse=", "), "\n")
cat("Unique regions:", paste(levels(main_data$Region), collapse=", "), "\n")
cat("Groups levels:", paste(levels(main_data$Groups), collapse=", "), "\n")

block_df <- main_data %>%
  group_by(subID, Region, Groups, property, distance) %>%
  summarise(mean_val = mean(OpticalProperty, na.rm = TRUE),
            n_obs = n(),
            .groups = "drop")

# Export basic counts for EDA
counts_table <- main_data %>%
  group_by(property, Region, Groups) %>%
  summarise(n_total = n(), n_blocks = n_distinct(subID), .groups = "drop")
write.csv(counts_table, file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", "counts_table.csv"), row.names = FALSE)

per_distance_results <- list()

props <- unique(block_df$property)
for(prop in props){
  dfp <- block_df %>% filter(property == prop)
  dists <- sort(unique(dfp$distance))
  res_list <- list()
  for(d in dists){
    subdf <- dfp %>% filter(distance == d)
    # collapse to subject mean across region if multiple regions per subject at same distance?
    # But we keep subID as unit
    # If only one observation per subID, lmer can still run but random effect may be redundant
    # Check unique subIDs per group
    nsub <- n_distinct(subdf$subID)
    # If fewer than 4 subjects, skip
    if(nsub < 4){
      res_list[[as.character(d)]] <- tibble(property = prop, distance = d,
                                            model = NA, p_group = NA, note = "too few subjects")
      next
    }
    # Fit LMM on mean_val
    # If model fails, fallback to t-test (unpaired)
    mod <- tryCatch(lmer(mean_val ~ Groups + (1 | subID), data = subdf, REML = TRUE), error = function(e) NULL)
    if(is.null(mod)){
      # fallback: t-test on group means
      tt <- tryCatch(t.test(mean_val ~ Groups, data = subdf), error = function(e) NULL)
      pv <- ifelse(is.null(tt), NA, tt$p.value)
      res_list[[as.character(d)]] <- tibble(property = prop, distance = d, model = "t-test", p_group = pv, note = "lmer_failed")
    } else {
      summ <- summary(mod)
      # attempt to extract p-value for Groups effect (lmerTest gives it)
      mod_anova <- anova(mod)
      pval <- if("Groups" %in% rownames(mod_anova)) mod_anova["Groups", "Pr(>F)"] else NA
      
      res_list[[as.character(d)]] <- tibble(property = prop, distance = d, model = "lmer", p_group = pval, note = NA)
    }
  }
  per_distance_results[[as.character(prop)]] <- bind_rows(res_list)
  write.csv(per_distance_results[[as.character(prop)]],
            file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", paste0("per_distance_results_", prop, ".csv")), row.names = FALSE)
}  #ERROR NEEDS FIXING

# 4. Global LMM (distance continuous) per property
#    mean_val ~ Groups * dist_c + (1|subID)
# -----------------------------
global_results <- list()

for(prop in props){
  datp <- block_df %>% filter(property == prop)
  datp <- datp %>% mutate(dist_c = distance - mean(distance, na.rm = TRUE))
  
  # Need enough subjects to estimate random effect
  if(n_distinct(datp$subID) < 4){
    global_results[[as.character(prop)]] <- tibble(property = prop, model = NA, note = "too few subjects")
    next
  }
  
  modg <- tryCatch(lmer(mean_val ~ Groups * dist_c + (1 | subID), data = datp, REML = TRUE),
                   error = function(e) NULL)
  
  if(is.null(modg)){
    global_results[[as.character(prop)]] <- tibble(property = prop, model = "lmer_failed", note = NA)
  } else {
    # Save summary and ANOVA to file
    sm <- summary(modg)
    an <- anova(modg)
    sink(file.path("~/Documents/MA 675/Consulting-Project-EPVs", paste0("global_lmm_summary_", prop, ".txt")))
    print(sm)
    print(an)
    sink()
    
    # Extract p-values from ANOVA
    p_group <- if("Groups" %in% rownames(an)) an["Groups", "Pr(>F)"] else NA
    p_inter <- if("Groups:dist_c" %in% rownames(an)) an["Groups:dist_c", "Pr(>F)"] else NA
    
    # Save results to the list
    global_results[[as.character(prop)]] <- tibble(
      property = prop,
      model = "lmer",
      p_group = p_group,
      p_interaction = p_inter,
      note = NA
    )
  }
}

# Combine into a single dataframe
global_results_df <- bind_rows(global_results)

# Save
write.csv(global_results_df, file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", "global_lmm_results.csv"), row.names = FALSE)


# 5. Diagnostics for global LMM (example: for scattering)
#    Save residuals vs fitted and QQ plots
# -----------------------------
diagnostic_plot_fn <- function(mod, label){
  # Resid vs Fitted
  dfres <- data.frame(fitted = fitted(mod), resid = resid(mod))
  p1 <- ggplot(dfres, aes(x = fitted, y = resid)) +
    geom_point() + geom_smooth(se = FALSE, method = "loess") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = paste0("Residuals vs Fitted: ", label))
  # QQ plot
  qq <- qqnorm(resid(mod), main = paste0("QQ plot: ", label)); qqline(resid(mod))
  # Save p1
  ggsave(file.path("~/Documents/MA 675/Consulting-Project-EPVs", paste0("resid_fitted_", label, ".png")), p1, width = 6, height = 4, dpi = 150)
  # For QQ, capture to png
  png(file.path("~/Documents/MA 675/Consulting-Project-EPVs", paste0("qqplot_", label, ".png")), width = 600, height = 400)
  qqnorm(resid(mod)); qqline(resid(mod))
  dev.off()
}
# Run diagnostics for each property model if available
for(prop in props){
  # try to read stored global model summary file to check if exists
  # Refit quickly to pass model object here (we already fitted earlier)
  datp <- block_df %>% filter(property == prop) %>% mutate(dist_c = distance - mean(distance, na.rm = TRUE))
  modg <- tryCatch(lmer(mean_val ~ Groups * dist_c + (1 | subID), data = datp, REML = TRUE), error = function(e) NULL)
  if(!is.null(modg)){
    diagnostic_plot_fn(modg, paste0("global_", prop))
  }
}

# 6. If linearity fails: GAMM (gamm4) with group-specific smooths
# -----------------------------
gamm_results <- list()
for(prop in props){
  datp <- block_df %>% filter(property == prop)
  if(nrow(datp) < 10){
    gamm_results[[as.character(prop)]] <- tibble(property = prop, note = "too few observations")
    next
  }
  # Fit GAMM with group-specific smooths s(distance, by=Groups)
  # Convert Groups to character for 'by'
  datp$Groups_char <- as.character(datp$Groups)
  gam_mod <- tryCatch(gamm4(mean_val ~ Groups_char + s(distance, by = Groups_char), random = ~(1|subID), data = datp), error = function(e) NULL)
  if(is.null(gam_mod)){
    gamm_results[[as.character(prop)]] <- tibble(property = prop, note = "gamm4_failed")
  } else {
    # Save summary
    sink(file.path("~/Documents/MA 675/Consulting-Project-EPVs", paste0("gamm4_summary_", prop, ".txt")))
    print(summary(gam_mod$gam))
    sink()
    gamm_results[[as.character(prop)]] <- tibble(property = prop, note = "gamm4_ok")
  }
}
write.csv(bind_rows(gamm_results), file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", "gamm4_status.csv"), row.names = FALSE)

# 7. Propensity score matching (subject-level)
#    We need covariates to match on. We'll use Region and subject-level baseline mean.
# -----------------------------
# Create subject-level baseline mean per property (aggregated)
subject_baseline <- main_data %>%
  group_by(subID, property, Region, Groups) %>%
  summarise(subject_mean = mean(OpticalProperty, na.rm = TRUE), .groups = "drop")
# ERROR ABOVE
psm_results <- list()
for(prop in props){
  sb <- subject_baseline %>% filter(property == prop)
  # ensure Groups is factor with two levels
  if(length(unique(sb$Groups)) < 2) {
    psm_results[[as.character(prop)]] <- tibble(property = prop, note = "not enough groups")
    next
  }
  # Fit propensity score model: predict Groups from Region + subject_mean
  # Convert Groups to binary for MatchIt
  sb_for_match <- sb %>% mutate(Groups = factor(as.character(Groups)))
  # if too few subjects, skip
  if(nrow(sb_for_match) < 6) {
    psm_results[[as.character(prop)]] <- tibble(property = prop, note = "too few subjects for match")
    next
  }
  m.out <- tryCatch(matchit(Groups ~ Region + subject_mean, data = sb_for_match, method = "nearest", ratio = 1), error = function(e) NULL)
  if(is.null(m.out)){
    psm_results[[as.character(prop)]] <- tibble(property = prop, note = "matchit_failed")
    next
  }
  matched <- match.data(m.out)
  # Save matched subject ids
  matched_ids <- matched$subID
  # Now use matched subjects to create block-level matched data
  matched_blocks <- block_df %>% filter(property == prop, subID %in% matched_ids)
  # Fit global LMM on matched_blocks
  mod_matched <- tryCatch(lmer(mean_val ~ Groups * (distance - mean(distance)) + (1 | subID), data = matched_blocks, REML = TRUE), error = function(e) NULL)
  if(is.null(mod_matched)){
    psm_results[[as.character(prop)]] <- tibble(property = prop, note = "model_failed_after_match")
  } else {
    # record p-values
    coefs <- coef(summary(mod_matched))
    p_group <- if("Groupsexperimental" %in% rownames(coefs)) coefs["Groupsexperimental","Pr(>|t|)"] else NA
    p_inter <- if("Groupsexperimental:(distance - mean(distance))" %in% rownames(coefs)) coefs["Groupsexperimental:(distance - mean(distance))","Pr(>|t|)"] else NA
    psm_results[[as.character(prop)]] <- tibble(property = prop, p_group_matched = p_group, p_inter_matched = p_inter, n_matched = length(unique(matched$subID)))
  }
}
psm_df <- bind_rows(psm_results)
write.csv(psm_df, file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", "psm_results.csv"), row.names = FALSE)

# 8. Balanced bootstrap to account for Group imbalance
#    We resample block-level observations, ensuring equal group sizes each iteration
# -----------------------------
bootstrap_results <- list()
B <- 500

for(prop in props){
  datp <- block_df %>% filter(property == prop)
  group_sizes <- table(datp$Groups)
  min_group_n <- min(group_sizes)
  
  if(min_group_n < 2){
    bootstrap_results[[as.character(prop)]] <- tibble(property = prop, note = "insufficient group size")
    next
  }
  
  coef_store <- numeric(B)
  p_store <- numeric(B)
  
  for(i in 1:B){
    samp <- datp %>%
      group_by(Groups) %>%
      sample_n(min_group_n, replace = TRUE) %>%
      ungroup()
    
    modb <- tryCatch(lmer(mean_val ~ Groups + (1 | subID), data = samp, REML = TRUE), error = function(e) NULL)
    
    if(is.null(modb)){
      coef_store[i] <- NA
      p_store[i] <- NA
    } else {
      # Safe extraction
      coefs <- coef(summary(modb))
      row_group <- grep("^Groups", rownames(coefs), value = TRUE)
      if(length(row_group) > 0 & "Pr(>|t|)" %in% colnames(coefs)){
        coef_store[i] <- coefs[row_group[1], "Estimate"]
        p_store[i] <- coefs[row_group[1], "Pr(>|t|)"]
      } else {
        coef_store[i] <- NA
        p_store[i] <- NA
      }
    }
  }
  
  bootstrap_results[[as.character(prop)]] <- tibble(
    property = prop,
    coef_median = median(coef_store, na.rm = TRUE),
    coef_mean = mean(coef_store, na.rm = TRUE),
    p_less_0.05 = mean(p_store < 0.05, na.rm = TRUE),
    n_boot = B
  )
}

bootstrap_df <- bind_rows(bootstrap_results)
write.csv(bootstrap_df, file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", "balanced_bootstrap_results.csv"), row.names = FALSE)


#9. Save summary outputs and make combined diagnostic plots
# -----------------------------
# Per-distance results and global results already saved.
# Also produce a visualization of per-distance p-values for each property
# Combine all per-distance results into one data frame
# Assuming per_distance_results is a list of dataframes (one per property)
all_per_dist <- bind_rows(per_distance_results)  # this will include property, distance, p_group, model, note

# Make sure p_group is numeric
all_per_dist$p_group <- as.numeric(all_per_dist$p_group)

# Save CSV
write.csv(all_per_dist,
          file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", "per_distance_all_properties.csv"),
          row.names = FALSE)

all_per_dist <- read.csv(file.path("~/Documents/MA 675/Consulting-Project-EPVs", "per_distance_all_properties.csv"))
ggplot(all_per_dist, aes(x = distance, y = p_group, color = property)) +
  geom_point() + geom_line() +
  facet_wrap(~property, scales = "free_y") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Per-distance p-values (lmer where possible)", y = "p-value") -> p_pvals
ggsave(file.path("~/Documents/MA 675/Consulting-Project-EPVs", "per_distance_pvalues.png"), p_pvals, width = 9, height = 4, dpi = 200)

# Save workspace for inspection
save.image(file = file.path("~/Documents/MA 675/Consulting-Project-EPVs", "analysis_workspace.RData"))

cat("All done. Outputs written to:", "~/Documents/MA 675/Consulting-Project-EPVs", "\n")
