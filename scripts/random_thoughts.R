test_requests <- tibble(period = c(2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024),
                        cases = c(NA, 1191, 1103, 927, 857, 946, 910, NA, NA),
                        tests = c(17485, 19008, 20025, 24702, 14184, 14761, 14754, 16806, 16482),
                        population = c(NA, 55619430, 55977178, 56230056, 56325961, 56554891, 57106398, NA, NA))

plot_grid(plotlist = list(
  test_requests %>%
  ggplot() +
  geom_line(aes(x = period, y = cases)) +
  theme_minimal() +
  labs(x = "Year",
       y = "Cases",
       title = "Cases by year"),
  test_requests %>%
  ggplot() +
  geom_line(aes(x = period, y = tests)) +
  theme_minimal() +
  labs(x = "Year",
       y = "Tests",
       title = "Tests by year"),
  test_requests %>%
    mutate(tpr = cases/tests) %>%
    ggplot() +
    geom_line(aes(x = period, y = tpr)) +
    theme_minimal() +
    labs(x = "Year",
         y = "Test-positivity rate (Cases/Samples tested)",
         title = "Test positivity rate by year")),
  ncol = 1)


# Tourism locations density -----------------------------------------------
pacman::p_load(osmdata, exactextractr, sf)

# Convert to OSM bounding box
aoi_osm <- opq("England") %>%
  add_osm_feature(key = "tourism", value = c("camp_pitch", "camp_site", "caravan_site", "picnic_site", "viewpoint")) %>%
  osmdata_sf()

osm_vect <- vect(st_transform(aoi_osm$osm_points, crs(gpw4_adjusted)))

campsite_raster <- rasterize(osm_vect, gpw4_adjusted, fun = "count", background = NA)

writeRaster(campsite_raster, here("data", "campsite_raster.tiff"))
