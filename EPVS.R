library(readxl)
library(dplyr)
library(ggplot2)
library(boot)

# Read in the dataset
scattering <- read_excel("~/Downloads/Datasets/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "scattering")
retardance <- read_excel("~/Downloads/Datasets/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "retardance")
orientation <- read_excel("~/Downloads/Datasets/caa_all_radii_40um_donut_13Oct2025.xlsx", sheet = "orientation")


scattering$property <- "scattering"
retardance$property <- "retardance"
orientation$property <- "orientation"

main_data <- bind_rows(scattering, retardance, orientation)



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


