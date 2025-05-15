###
# ONS population data
###

pacman::p_load(onsr, tidyverse, disaggregation, terra, tidyterra, here, cowplot)

# Lyme data
data <- read_csv(here("data", "FULL-UKHSA-2017-2022-Lyme-Disease.csv"))

lad_lyme <- unique(data$`Area Code`)

# around 1km resolution
gpw4 <- rast(here("large-data", "GPW4.tiff"))

# Council mapping
lad <- vect(here("large-data", "LAD_DEC_24_UK_BFC.shp")) %>%
  simplifyGeom(tolerance = 20)

writeVector(lad, here("data", "reduced_LAD_boundaries.rds"))

# County lyme data
data_subset <- data %>%
  filter(`Area Code` %in% lad$LAD24CD)

data_summary <- data_subset %>%
  group_by(`Area Code`, `Time period`) %>%
  summarise(count = sum(Count, na.rm = TRUE),
            population = sum(Denominator, na.rm = TRUE),
            period = `Time period`) %>%
  ungroup()

data_wide <- data_summary %>%
  mutate(
    count_col = paste0("count_", period),
    pop_col = paste0("population_", period)
  ) %>%
  # Pivot wider twice, once for count and once for population
  select(`Area Code`, period, count, population) %>%
  pivot_wider(
    names_from = period,
    values_from = c(count, population),
    names_glue = "{.value}_{period}"
  )

lad_with_cases <- lad %>%
  left_join(data_wide, by = c("LAD24CD" = "Area Code"))

writeVector(lad_with_cases, here("data", "LAD_boundaries_with_cases.rds"))

lad_joined <- left_join(lad, data_summary, by = c("LAD24CD" = "Area Code"), copy = TRUE)

lad_plots <- list()

periods <- unique(lad_joined$period)

for (p in periods[1:6]) {
  
  lad_subset <- lad_joined[lad_joined$period == p, ]
  
  p_plot <- ggplot(lad_subset) +
    geom_spatvector(aes(fill = count)) +
    scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
    labs(title = paste("Period:", p)) +
    theme_minimal()
  
  lad_plots[[as.character(p)]] <- p_plot
}

plot_grid(plotlist = lad_plots)

# Higher resolution pop raster match to the county polygons
lad_extent <- ext(lad)

lad_extent_gpw4_crs <- project(lad, crs(gpw4))

gpw4_crop <- crop(gpw4, lad_extent_gpw4_crs)

gpw4_crop_proj <- project(gpw4_crop, crs(lad))

lad_pop <- extract(gpw4_crop_proj, lad, fun = sum, na.rm = TRUE) %>%
  mutate(LAD24CD = lad$LAD24CD)

compare_to_data <- data_subset %>%
  left_join(lad_pop, by = c(`Area Code` = "LAD24CD")) %>%
  ggplot(aes(x = Denominator, y = GPW4)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +  
  coord_fixed(ratio = 1, xlim = c(0, 2000000), ylim = c(0, 2000000)) +          
  labs(x = "Data Denominator", y = "GPW4 Population Estimate",
       title = "Comparison of LAD Population Totals: GPW4 vs Provided Data") +
  theme_minimal()

# regression based adjustment of GPW4 data
lad_compare <- data_subset %>%
  left_join(lad_pop, by = c("Area Code" = "LAD24CD")) %>%
  rename(obs_pop = Denominator, gpw4_pop = GPW4)

lm_fit <- lm(obs_pop ~ gpw4_pop, data = lad_compare)
summary(lm_fit)

beta0 <- coef(lm_fit)[1]
beta1 <- coef(lm_fit)[2]

gpw4_corrected <- beta0 + beta1 * gpw4_crop_proj
names(gpw4_corrected) <- "gpw4_corrected"

lad_pop_corrected <- extract(gpw4_corrected, lad, fun = sum, na.rm = TRUE) %>%
  mutate(LAD24CD = lad$LAD24CD)

lad_compare <- lad_compare %>%
  select(-ID) %>%
  left_join(lad_pop_corrected, by = c("Area Code" = "LAD24CD"))

lm_fit <- lm(obs_pop ~ gpw4_corrected, data = lad_compare)
summary(lm_fit)

# Plot updated comparison
ggplot(lad_compare, aes(x = obs_pop, y = gpw4_corrected)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  coord_fixed(ratio = 1, xlim = c(0, 2000000), ylim = c(0, 2000000)) +
  labs(x = "Census Population", y = "Corrected GPW4 Population",
       title = "Post-Adjustment Comparison: GPW4 vs Census") +
  theme_minimal()
