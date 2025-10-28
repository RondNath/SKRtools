#' @title SKR analysis - compute regression parameters of the SKR
#' @description Launch the SKR analysis of datasets with skewness and kurtosis values, based on explanatory factors (e.g., spatial, temporal, context, ...)
#' @param moments dataframe of distribution moments with skewness and kurtosis columns
#' @param Factors dataframe of factors for the computation of statistics
#' @param lin_mod Indicates the type of linear model to use for (SKR): choose "lm" or "mblm"
#' @param slope_ref slope of a reference SKR used as a baseline (default: slope_ref = 1; skew-uniform slope)
#' @param intercept_ref intercept of a reference SKR used as a baseline (default: intercept_ref = 1.86; skew-uniform intercept)
#' @param distance_metric Indicates the method to compute distance-based regression parameters: choose "RMSE" (for Root Mean Square Error, default) or "MAE" (for Mean Absolute Error)
#' @param workers Cores used for parallelisation
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
#' @importFrom dplyr mutate select filter across slice %>%
#' @importFrom future plan
#' @importFrom future.apply future_lapply
#' @importFrom progressr with_progress progressor handlers
#' @importFrom mblm mblm
#' @importFrom parallel detectCores
#' @importFrom stats lm residuals var
#' @importFrom tidyr unite
#' @examples
#' head(SK)
#' SKRparam(
#' moments = SK,
#' Factors = SK[,c("Factor")],
#' slope_ref = 1,
#' intercept_ref = 1.86,
#' distance_metric = "RMSE",
#' lin_mod = "lm",
#' workers = parallel::detectCores() - 1
#' )

SKRparam <- function(
    moments,
    Factors = NULL,
    lin_mod = "lm",
    slope_ref = 1,
    intercept_ref = 1.86,
    distance_metric = "RMSE",
    workers = parallel::detectCores() - 1
) {
  # test data input validity
  if (is.null(moments) || is.null(moments$skewness) || is.null(moments$kurtosis)){
    stop("We need skewness and kurtosis values in 'moments'!")
  }
  if (lin_mod != "lm" && lin_mod != "mblm"){
    stop("Unknown lin_mod. Use 'lm' or 'mblm'!")
  }
  if (distance_metric != "RMSE" && distance_metric != "MAE"){
    stop("Unknown distance_metric. Use 'RMSE' or 'MAE'!")
  }
  if (!is.null(Factors) && (nrow(data.frame(Factors)) != nrow(data.frame(moments)))){
    stop("'Factors' and 'moments' must have the same lenght!")
  }
  if (!is.null(Factors) && any((data.frame(Factors) %>% 
          tidyr::unite("Factors", 
                       dplyr::everything(), 
                       sep = "_", 
                       remove = TRUE) %>%
          dplyr::count(Factors))$n < 2)) {
    stop("Must have at least two repetitions per Factors to build SKR!")
  }
  if (!is.null(Factors) && any((data.frame(Factors) %>% 
                                tidyr::unite("Factors", 
                                             dplyr::everything(), 
                                             sep = "_", 
                                             remove = TRUE) %>%
                                dplyr::count(Factors))$n == 2)) {
    warning("Warning: for some Factors, SKR is built with only two observations")
  }
  
  # session type for parallelisation
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, 
                 workers = workers)
  } else {
    future::plan(future::multicore, 
                 workers = workers)
  }
  
  # group Factors to build SKR
  if (!is.null(Factors)){
    moments <- cbind(moments,
                     data.frame(Factors) %>% 
                       tidyr::unite("Factors", 
                                    dplyr::everything(), 
                                    sep = "_", 
                                    remove = TRUE))
  }else{
    moments <- moments %>% 
      dplyr::mutate(Factors = "A")
  }
  return_SKR <- unique(moments %>%
                         dplyr::select(Factors, 
                                       colnames(Factors)))
  group_factor <- unique(moments$Factors)
  
  # progression bar
  results <- NULL
  progressr::with_progress({
    handlers("txtprogressbar")
    p <- progressr::progressor(along = group_factor)
    res_list <- future.apply::future_lapply(group_factor, 
                              function(i) {
      p()# progression
      
                                moments_filter <- moments %>% 
        dplyr::filter(Factors == i)
      y <- moments_filter$kurtosis
      x <- moments_filter$skewness^2
      distance_reference_SKR <- y - (slope_ref * x + intercept_ref)
      # linear model
      if (lin_mod == "lm") {
        fit <- lm(y ~ x)
      } else if (lin_mod == "mblm") {
        fit <- mblm::mblm(y ~ x)
      }
      # Compute SKR parameters
      if (distance_metric == "RMSE") {
        dist_pred <- sqrt(mean(fit$residuals^2, na.rm = TRUE))
        dist_ref <- sqrt(mean(distance_reference_SKR^2, na.rm = TRUE))
      } else if (distance_metric == "MAE") {
        dist_pred <- mean(abs(fit$residuals), na.rm = TRUE)
        dist_ref <- mean(abs(distance_reference_SKR), na.rm = TRUE)
      } else {
        stop("Unknown distance_metric. Use 'RMSE' or 'MAE'.")
      }
      cv_pred <- sd(abs(fit$residuals), na.rm = TRUE) * 100 / mean(abs(fit$residuals), na.rm = TRUE)
      cv_ref <- sd(abs(distance_reference_SKR), na.rm = TRUE) * 100 / mean(abs(distance_reference_SKR), na.rm = TRUE)
      
      # compil SKR parameters per Factors
      data.frame(
        return_SKR %>% 
          dplyr::filter(Factors == i) %>% 
          slice(1),
        Slope = coef(fit)[2],
        Intercept = coef(fit)[1],
        Rsquare = 1 - (mean(stats::residuals(fit)^2, na.rm = TRUE) / var(y, na.rm = TRUE)),
        distance_predicted_SKR = dist_pred,
        distance_reference_SKR = dist_ref,
        CV_distance_predicted_SKR = cv_pred,
        CV_distance_reference_SKR = cv_ref,
        stringsAsFactors = FALSE
      )
    })
    results <- do.call(rbind, res_list)
  })
  return(results)
}
