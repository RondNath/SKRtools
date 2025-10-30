# I. PACKAGES ----
library(devtools)
devtools::install_github("RondNath/SKRtools", force = T)
library(SKRtools)
library(dplyr)

# II. Compute moments ----
## Case 1: Compute moment for a list of values ("values") ----
SKRtools::Moments(
  values = SKRtools::values$Values,
  individual_id = NULL,
  weight = NULL,
  values_factors = NULL,
  weight_factors = NULL
)

# Case 2: Compute moments of a list of values ("values") observed under different conditions ("values_factors") ----
SKRtools::Moments(
  values = SKRtools::values$Values,
  individual_id = NULL,
  weight = NULL,
  values_factors = data.frame(SKRtools::values[,c("Sub_Factor", "Factor"), drop = FALSE]),
  weight_factors = NULL
)

SKRtools::Moments(
  values = SKRtools::values$Values,
  individual_id = NULL,
  weight = NULL,
  values_factors = data.frame(SKRtools::values[,c("Factor"), drop = FALSE]),
  weight_factors = NULL
)

# Case 3: Compute weighting ("weight" & "weight_factors") moments with or without variability in the observation of values by individuals ("values_factors") ---- 
SKRtools::Moments(
  values = (SKRtools::values %>% 
    filter(Factor == "A" & Sub_Factor == 1))$Values,
  individual_id = (SKRtools::values %>% 
                     filter(Factor == "A" & Sub_Factor == 1))$Individuals,
  weight = SKRtools::weight %>% 
    dplyr::select(-Factor, -Sub_Factor),
  values_factors = NULL,
  weight_factors = data.frame(SKRtools::weight[,c("Factor", "Sub_Factor")])
)

SKRtools::Moments(
  values = SKRtools::values$Values,
  individual_id = SKRtools::values$Individuals,
  weight = SKRtools::weight %>% 
    dplyr::select(-Factor, -Sub_Factor),
  values_factors = data.frame(SKRtools::values[,c("Factor", "Sub_Factor")]),
  weight_factors = data.frame(SKRtools::weight[,c("Factor", "Sub_Factor")])
)

# III. SKR analysis ----
SKRtools::SKRparam(
    moments = SKRtools::SK,
    Factors = SKRtools::SK[,c("Sub_Factor", "Factor"), drop = FALSE],
    lin_mod = "lm",
    slope_ref = 1,
    intercept_ref = 1.86,
    distance_metric = "RMSE"
)

SKRtools::SKRparam(
  moments = SKRtools::SK,
  Factors = NULL,
  lin_mod = "lm",
  slope_ref = 1,
  intercept_ref = 1.86,
  distance_metric = "RMSE"
)
