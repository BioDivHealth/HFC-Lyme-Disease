#---------------------------#
# Disaggregation regression #
#---------------------------#

# Load packages ----
pacman::p_load(disaggregation, tidyverse)

# Load & wrangle data ----

## Load shape file for LD incidence data aggregated to county level
# cols: 
# geometry
# id_var: name of column with polygon id
# response_var: response count data
agg_dat <- sf::read_sf("data/XXX.shp") 

## SpatRaster of covariate rasters to be used in the model
# One layer per covariate, at same spatial resolution for each covariate
covar_stack <- terra::SpatRaster("path/to/covariate_raster.tif")

## Define population raster at the same spatial scale as covariates
# Ensures correct aggregation from pixel -> polygon
# Not supplying this assumes uniform population across polygon
pop_raster <- XXX

# Prepare data for input into model
model_data <- prepare_data(
  polygon_shapefile = agg_dat, 
  covariate_rasters = covar_stack,
  aggregation_raster = pop_raster,
  id_var = "XXX",
  response_var = "XXX",
  # parameters for mesh to create spatial field
  mesh_args = list(cutoff = 0.01,
                   offset = c(0.1, 0.5),
                   max.edge = c(0.1, 0.2),
                   resolution = 250),
  # NAs in response or covariates give errors
  na_action = FALSE
  ...)

# Check what the data look like
plot(model_data)

# Fit the model & estimate parameters @ county level ----
model <- disag_model(
  data = model_data,   # prepared input data
  family = "poisson",  # likelihood function
  link = "log",        # link function
  iterations = 100,    # no. iterations for optimisation
  # named list of prior values (can tweak)
  priors = list(priormean_intercept = 0,
                priorsd_intercept = 2,
                priormean_slope = 0.0,
                priorsd_slope = 0.4,
                prior_rho_min = 3,
                prior_rho_prob = 0.01,
                prior_sigma_max = 1,
                prior_sigma_prob = 0.01,
                prior_iideffect_sd_max = 0.05,
                prior_iideffect_sd_prob = 0.01)
  )

# Make predictions for incidence at fine scale using fitted parameters ----

# Predict incidence at scale of covariates
preds <- predict(model, N = 100, CI = 0.95)

# Visualise fine-scale incidence predictions
plot(preds)

# Some step to aggregate these fine scale predictions back up to county level
## ~~~ code here ~~~~

