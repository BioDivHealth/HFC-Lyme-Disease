
# Load packages -----------------------------------------------------------
pacman::p_load(onsr, tidyverse, disaggregation, terra, tidyterra, here, cowplot)


# Load data ---------------------------------------------------------------

# Lyme data
data <- read_csv(here("data", "FULL-UKHSA-2017-2022-Lyme-Disease.csv"))

lad_lyme <- unique(data$`Area Code`)

# around 1km resolution
gpw4 <- rast(here("large-data", "GPW4.tiff"))

# Council mapping
lad <- vect(here("large-data", "LAD_DEC_24_UK_BFC.shp")) %>%
  simplifyGeom(tolerance = 20)

writeVector(lad, here("data", "reduced_LAD_boundaries.rds"))


# Process data for summaries ----------------------------------------------

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


# Correct the GPW4 data to more closely match observed populations --------

lad_extent <- ext(lad)

lad_extent_gpw4_crs <- project(lad, crs(gpw4))

gpw4_crop <- crop(gpw4, lad_extent_gpw4_crs)

gpw4_crop_proj <- project(gpw4_crop, crs(lad))

lad_pop <- extract(gpw4_crop_proj, lad, fun = sum, na.rm = TRUE, bind = TRUE) %>%
  as_tibble() %>%
  select(LAD24CD, GPW4) 

lad_compare <- data_subset %>%
  rename(obs_pop = Denominator) %>%
  left_join(lad_pop, by = c("Area Code" = "LAD24CD")) %>%
  mutate(scaling_factor = obs_pop / GPW4)

lad_sf <- lad %>%
  left_join(lad_compare %>% select(LAD24CD = `Area Code`, scaling_factor), by = "LAD24CD")

scaling_raster <- rasterize(lad_sf, gpw4_crop_proj, field = "scaling_factor")

gpw4_adjusted <- gpw4_crop_proj * scaling_raster

plot(crop(gpw4_adjusted, lad))

# Compare GPW4 now to observed data
lad_pop_adjusted <- extract(gpw4_adjusted, lad, fun = sum, na.rm = TRUE, bind = TRUE) %>%
  as_tibble() %>%
  select(LAD24CD, GPW4)

lad_compare_models <- data_subset %>%
  select(`Area Code`, Denominator) %>%
  rename(obs_pop = Denominator) %>%
  left_join(lad_pop %>%
              rename(GPW4_original = GPW4), by = c("Area Code" = "LAD24CD")) %>%
  left_join(lad_pop_adjusted %>%
              rename(GPW4_scaled = GPW4), by = c("Area Code" = "LAD24CD"))

model_original <- lm(obs_pop ~ GPW4_original, data = lad_compare_models)
model_adjusted <- lm(obs_pop ~ GPW4_scaled, data = lad_compare_models)

r2_original <- summary(model_original)$r.squared
r2_adjusted <- summary(model_adjusted)$r.squared

plot_original <- ggplot(lad_compare_models, aes(x = GPW4_original, y = obs_pop)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "Original GPW4 vs Observed",
       subtitle = paste0("R² = ", round(r2_original, 3)),
       x = "GPW4 Population (Original)",
       y = "Observed Population") +
  coord_fixed() +
  theme_minimal()

plot_adjusted <- ggplot(lad_compare_models, aes(x = GPW4_scaled, y = obs_pop)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  labs(title = "Adjusted GPW4 vs Observed",
       subtitle = paste0("R² = ", round(r2_adjusted, 3)),
       x = "GPW4 Population (Adjusted)",
       y = "Observed Population") +
  coord_fixed() +
  theme_minimal()

plot_grid(plotlist = list(plot_original, plot_adjusted))

writeRaster(gpw4_adjusted, here("data", "gpw4_adjusted_raster.tiff"))
