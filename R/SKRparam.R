#' @title SKR analysis - compute regression parameters of the SKR
#' @description Launch the SKR analysis of datasets with skewness and kurtosis values, based on explanatory factors (e.g., spatial, temporal, context, ...)
#' @param moments dataframe of distribution moments with skewness and kurtosis columns
#' @param Factors dataframe of factors for the computation of statistics
#' @param lin_mod Indicates the type of linear model to use for (SKR): choose "lm" or "mblm"
#' @param slope_ref slope of a reference SKR used as a baseline (default: slope_ref = 1; skew-uniform slope)
#' @param intercept_ref intercept of a reference SKR used as a baseline (default: intercept_ref = 1.86; skew-uniform intercept)
#' @param distance_metric Indicates the method to compute distance-based regression parameters: choose "RMSE" (for Root Mean Square Error, default) or "MAE" (for Mean Absolute Error)
#' @returns
#' data.frame with:
#' Slope per Factors,
#' Intercept per Factors,
#' Rsquare per Factors,
#' distance from predicted SKR per Factors,
#' CV of the distance from predicted SKR per Factors.
#' distance from reference SKR per Factors,
#' CV of the distance from reference SKR per Factors.
#' @export
#' @importFrom data.table as.data.table .SD := .N
#' @importFrom stats lm residuals var
#' @importFrom mblm mblm
#' @examples
#' head(SK)
#' SKRparam(
#' moments = SK,
#' Factors = SK[,c("Factor")],
#' slope_ref = 1,
#' intercept_ref = 1.86,
#' distance_metric = "RMSE",
#' lin_mod = "lm"
#' )

SKRparam <- function(
    moments,
    Factors = NULL,
    lin_mod = "lm",
    slope_ref = 1,
    intercept_ref = 1.86,
    distance_metric = "RMSE"
) {
  # test data input validity
  if (is.null(moments) || is.null(moments$skewness) || is.null(moments$kurtosis)) {
    stop("We need skewness and kurtosis values in 'moments'!")
  }
  if (!lin_mod %in% c("lm", "mblm")) stop("Unknown lin_mod. Use 'lm' or 'mblm'!")
  if (!distance_metric %in% c("RMSE", "MAE")) stop("Unknown distance_metric. Use 'RMSE' or 'MAE'!")
  
  # intern function to do SKR analysis and compute parameters
  skr_analysis <- function(dt) {
    y <- dt$kurtosis
    x <- dt$skewness^2
    dist_ref <- y - (slope_ref * x + intercept_ref)
    
    fit <- if (lin_mod == "lm") lm(y ~ x) else mblm::mblm(y ~ x)
    
    residuals_fit <- stats::residuals(fit)
    dist_pred <- if (distance_metric == "RMSE") sqrt(mean(residuals_fit^2, na.rm = TRUE)) else mean(abs(residuals_fit), na.rm = TRUE)
    dist_ref_val <- if (distance_metric == "RMSE") sqrt(mean(dist_ref^2, na.rm = TRUE)) else mean(abs(dist_ref), na.rm = TRUE)
    
    cv_pred <- stats::sd(abs(residuals_fit), na.rm = TRUE) * 100 / mean(abs(residuals_fit), na.rm = TRUE)
    cv_ref <- stats::sd(abs(dist_ref), na.rm = TRUE) * 100 / mean(abs(dist_ref), na.rm = TRUE)
    
    data.table(
      Slope = stats::coef(fit)[2],
      Intercept = stats::coef(fit)[1],
      Rsquare = 1 - mean(residuals_fit^2, na.rm = TRUE) / stats::var(y, na.rm = TRUE),
      distance_predicted_SKR = dist_pred,
      distance_reference_SKR = dist_ref_val,
      CV_distance_predicted_SKR = cv_pred,
      CV_distance_reference_SKR = cv_ref
    )
  }
  
  dt_moments <- data.table::as.data.table(moments)
  
  if (!is.null(Factors)) {
    if ((nrow(Factors) != nrow(moments))) stop("'Factors' and 'moments' must have the same lenght!")
    groups <- colnames(Factors)
    if (any(dt_moments[, .N, by = groups]$N < 2)) {
    stop("Must have at least two observations per Factors to build SKR!")
    }
    res_SKR <- dt_moments[, skr_analysis(.SD), by = groups]
  } else {
    res_SKR <- skr_analysis(dt_moments)
  }
  return(res_SKR)
}