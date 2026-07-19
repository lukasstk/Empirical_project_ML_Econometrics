# =============================================================================
# 01_data_preparation.R
# Load the panel and construct outcomes / treatments / control sets.
# =============================================================================

source("Code/00_setup.R")

# load() restores the data frame under the name it was saved with: `dataset`
load("Data/urban_emissions_panel.RData")
data <- dataset
rm(dataset)

# ---- (1) Outcome: log emissions ----------------------------------------------
data$log_transport_co2 <- log(data$transport_co2)

# ---- (2) Treatments and their interaction ------------------------------------
data$cp_x_lez <- data$cp_active * data$lez_active

# ---- (3) Transformed controls -------------------------------------------------
# Logs of skewed level variables, factor coding for categorical variables
data$log_population  <- log(data$population)
data$log_pop_density <- log(data$pop_density)
data$log_gdp_pc      <- log(data$gdp_pc)
data$log_area_km2    <- log(data$area_km2)
data$log_elevation         <- log(data$elevation)
data$log_tourism_intensity <- log(data$tourism_intensity)

data$latitude_zone <- factor(data$latitude_zone)
data$country_id    <- factor(data$country_id)

# ---- (4) Control sets ----------------------------------------------------------
# Baseline controls = all columns minus the exclusion list below.
all_cols <- colnames(data)

exclude_from_controls <- c(
  # panel identifiers / structure (year and country_id enter via build_W)
  "city_id", "year", "country_id",

  # outcome and its raw pre-log form
  "transport_co2", "log_transport_co2",

  # treatment indicators and their timing variables
  "cp_active", "lez_active", "cp_x_lez",
  "cp_impl_year", "lez_impl_year", "cp_announce_year", "lez_announce_year",

  # mediators / bad controls (re-added as ctrl_mediators in the sensitivity
  # check, except total_co2 which contains the outcome)
  "total_co2", "pm25",
  "fleet_diesel_share", "fleet_petrol_share", "fleet_electric_share",

  # omitted reference category (industry shares sum to 1)
  "industry_public",

  # raw levels superseded by their log versions
  "population", "pop_density", "gdp_pc", "area_km2",
  "elevation", "tourism_intensity"
)

ctrl_baseline <- setdiff(all_cols, exclude_from_controls)

# Mediators for the sensitivity check (fleet_petrol_share dropped: shares sum to 1)
ctrl_mediators <- c("pm25", "fleet_diesel_share", "fleet_electric_share")

# Squared terms for key continuous controls (added via I(x^2) in build_W)
sq_vars <- c("log_population", "log_gdp_pc", "fuel_price",
             "public_transit_score")

# ---- (5) Save prepared objects for the analysis scripts ---------------------
saveRDS(list(data = data,
             ctrl_baseline  = ctrl_baseline,
             ctrl_mediators = ctrl_mediators,
             sq_vars        = sq_vars),
        file.path(out_dir, "prepared_data.rds"))

cat("\nData preparation done:", nrow(data), "city-years,",
    length(unique(data$city_id)), "cities.\n")
