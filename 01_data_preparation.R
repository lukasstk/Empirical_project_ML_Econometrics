# =============================================================================
# 01_data_preparation.R
# Load the panel, construct outcomes / treatments / control sets, and
# produce descriptive figures and tables.
# Corresponds to report section "Data preparation".
# =============================================================================

source("00_setup.R")

# load() restores the data frame under the name it was saved with: `dataset`
load("urban_emissions_panel.RData")
data <- dataset
rm(dataset)

# ---- (1) Outcomes -----------------------------------------------------------
# Log transformation: effects are interpreted as (approximate) percentage
# changes in emissions, and the strongly right-skewed emissions distribution
# becomes better behaved for the ML learners.
data$log_transport_co2 <- log(data$transport_co2)

# ---- (2) Treatments ---------------------------------------------------------
# cp_active and lez_active are the policy indicators.
# The product term allows the joint effect to differ from the sum of the
# individual effects (question 3: super-/sub-additivity).
data$cp_x_lez <- data$cp_active * data$lez_active

# ---- (3) Transformed controls ----------------------------------------------
# Logs of strongly skewed level variables
data$log_population  <- log(data$population)
data$log_pop_density <- log(data$pop_density)
data$log_gdp_pc      <- log(data$gdp_pc)
data$log_area_km2    <- log(data$area_km2)

data$latitude_zone <- factor(data$latitude_zone)
data$country_id    <- factor(data$country_id)

# ---- (4) Control sets -------------------------------------------------------
# BASELINE controls: potential confounders that can drive both policy
# adoption and emissions (size, economy, politics, prices, weather,
# geography, industry structure).
#
# Deliberately ALSO included: variables that should be irrelevant
# (museum visitors, benches, flagpoles, fountains, sister cities, ...).
# The lasso should set their coefficients to zero - a built-in check that
# variable selection works.
#
# Deliberately EXCLUDED from the baseline (bad controls):
#   - pm25:                air pollution is itself an outcome of emissions
#   - fleet_*_share:       fleet composition is the main MECHANISM through
#                          which a LEZ works (mediator, not confounder)
#   - total_co2:           contains the outcome variable
# They are added back in the sensitivity analysis to show what happens.
#
# Note: industry shares sum to 1, so industry_public is dropped to avoid
# perfect collinearity (it is the omitted reference category).
ctrl_baseline <- c(
  "log_population", "log_pop_density", "log_gdp_pc", "unemployment",
  "education_share", "public_transit_score", "road_km_pc",
  "tourism_intensity", "logistics_activity",
  "political_green", "fiscal_capacity", "electoral_competitiveness",
  "ngo_environment_index", "national_climate_pact",
  "electricity_price", "fuel_price",
  "temp_anomaly", "precip", "heating_degree_days",
  "renewable_electricity_share",
  "industry_manufacturing", "industry_services", "industry_logistics",
  "log_area_km2", "elevation", "coastal", "latitude_zone",
  # placebo / noise variables (lasso should discard these):
  "museum_visitors_pc", "library_count", "streetlight_density",
  "fountain_count", "bench_count_pc", "flagpole_count", "sister_city_count"
)

# Mediators / bad controls for the sensitivity check.
# fleet_petrol_share dropped: the three fleet shares sum to 1.
ctrl_mediators <- c("pm25", "fleet_diesel_share", "fleet_electric_share")

# Squared terms for key continuous controls (added via I(x^2) in build_W)
sq_vars <- c("log_population", "log_gdp_pc", "fuel_price",
             "public_transit_score")

# ---- (5) Heterogeneity variables (question 2) -------------------------------
# City characteristics along which the policy effects may differ.
# They are CENTERED before interacting with the treatments, so that the
# main treatment coefficient is the effect for an "average" city
# (this makes the baseline coefficient directly interpretable - unlike in
# the wage-gap example, where the main coefficient was a reference-group gap).
het_vars <- c("log_population", "log_gdp_pc", "log_pop_density",
              "public_transit_score", "political_green", "fuel_price",
              "tourism_intensity", "industry_logistics", "coastal")
for (v in het_vars) data[[paste0(v, "_c")]] <- data[[v]] - mean(data[[v]])

# ---- (6) Descriptives -------------------------------------------------------

# 6a. Policy adoption over time
adoption <- data %>%
  group_by(year) %>%
  summarise(`Congestion pricing` = sum(cp_active),
            `Low-emission zone`  = sum(lez_active),
            `Both policies`      = sum(cp_active * lez_active),
            .groups = "drop") %>%
  pivot_longer(-year, names_to = "policy", values_to = "n_cities")

p_adopt <- ggplot(adoption, aes(year, n_cities, colour = policy)) +
  geom_line(linewidth = 1) + geom_point(size = 1.5) +
  labs(x = "Year", y = "Number of cities with active policy",
       title = "Staggered policy adoption, 2008-2024", colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("output/figures/fig_adoption.png", p_adopt, width = 8, height = 5)

# 6b. Raw emission trends by (eventual) treatment group
data <- data %>%
  group_by(city_id) %>%
  mutate(group = case_when(
    max(cp_active) == 1 & max(lez_active) == 1 ~ "Both (eventually)",
    max(cp_active) == 1                        ~ "CP only",
    max(lez_active) == 1                       ~ "LEZ only",
    TRUE                                       ~ "Never treated")) %>%
  ungroup()

trends <- data %>%
  group_by(year, group) %>%
  summarise(mean_log_co2 = mean(log_transport_co2), .groups = "drop")

p_trend <- ggplot(trends, aes(year, mean_log_co2, colour = group)) +
  geom_line(linewidth = 1) +
  labs(x = "Year", y = "Mean log transport CO2 (kt)",
       title = "Raw emission trends by eventual policy status",
       colour = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "bottom")
ggsave("output/figures/fig_raw_trends.png", p_trend, width = 8, height = 5)

# 6c. Summary statistics of the key variables
sum_vars <- c("transport_co2", "cp_active", "lez_active",
              "population", "gdp_pc", "pop_density",
              "public_transit_score", "fuel_price", "political_green")
sumstats <- data.frame(
  variable = sum_vars,
  mean = sapply(sum_vars, function(v) mean(data[[v]])),
  sd   = sapply(sum_vars, function(v) sd(data[[v]])),
  min  = sapply(sum_vars, function(v) min(data[[v]])),
  max  = sapply(sum_vars, function(v) max(data[[v]]))
)
save_table(sumstats, "tab_summary_stats")

# 6d. Treatment counts
counts <- data.frame(
  quantity = c("cities total", "cities ever CP", "cities ever LEZ",
               "cities ever both", "city-years CP active",
               "city-years LEZ active", "city-years both active"),
  value = c(length(unique(data$city_id)),
            length(unique(data$city_id[data$cp_active == 1])),
            length(unique(data$city_id[data$lez_active == 1])),
            length(unique(data$city_id[data$cp_active == 1 & data$lez_active == 1])),
            sum(data$cp_active), sum(data$lez_active), sum(data$cp_x_lez))
)
save_table(counts, "tab_treatment_counts")

# ---- (7) Save prepared objects for the analysis scripts ---------------------
saveRDS(list(data = data,
             ctrl_baseline  = ctrl_baseline,
             ctrl_mediators = ctrl_mediators,
             sq_vars        = sq_vars,
             het_vars       = het_vars),
        "output/prepared_data.rds")

cat("\nData preparation done:", nrow(data), "city-years,",
    length(unique(data$city_id)), "cities.\n")
