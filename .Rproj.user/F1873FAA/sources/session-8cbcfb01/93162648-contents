# I. PACKAGES ----
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
  weight_factors = NULL,
  workers = parallel::detectCores() - 1
)

# Case 2: Compute moments of a list of values ("values") observed under different conditions ("values_factors") ----
SKRtools::Moments(
  values = SKRtools::values$Values,
  individual_id = NULL,
  weight = NULL,
  values_factors = SKRtools::values[,c("Sub_Factor"), drop = FALSE],
  weight_factors = NULL,
  workers = parallel::detectCores() - 1
)

# Case 3: Compute weighting ("weight" & "weight_factors") moments with or without variability in the observation of values by individuals ("values_factors") ---- 
SKRtools::Moments(
  values = SKRtools::values$Values,
  individual_id = SKRtools::values$Individuals,
  weight = SKRtools::weight %>% 
    dplyr::select(-Factor, -Sub_Factor),
  values_factors = NULL,
  weight_factors = SKRtools::weight[,c("Factor", "Sub_Factor")],
  workers = parallel::detectCores() - 1
)

SKRtools::Moments(
  values = SKRtools::values$Values,
  individual_id = SKRtools::values$Individuals,
  weight = SKRtools::weight %>% 
    dplyr::select(-Factor, -Sub_Factor),
  values_factors = SKRtools::values[,c("Factor", "Sub_Factor")],
  weight_factors = SKRtools::weight[,c("Factor", "Sub_Factor")],
  workers = parallel::detectCores() - 1
)

# III. SKR analysis ----
SKRtools::SKRparam(
    moments = SKRtools::SK,
    Factors = SKRtools::SK[,"Sub_Factor", drop = FALSE],
    lin_mod = "lm",
    slope_ref = 1,
    intercept_ref = 1.86,
    distance_metric = "RMSE",
    workers = parallel::detectCores() - 1
)