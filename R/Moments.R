#' @title SKR analysis - compute regression parameters of the SKR
#' @description Launch the SKR analysis of datasets with skewness and kurtosis values, based on explanatory factors (e.g., spatial, temporal, context, ...)
#' @param values vector of values observed/measured per individual (length format)
#' @param individual_id vector of Individual name (length format, same length as values, must be completed if weight =! NULL and correspond to the column name in the weighting data frame)
#' @param weight dataframe of the weight to attribute to each individual for each factors of observation (e.g., frequency, abundance, ..., NO NA) (wide format, column name = Individual observed, must correspond to the row name in the values data frame)
#' @param values_factors dataframe of the factors of observations of the values vector
#' @param weight_factors dataframe of the factors of observations of the weight dataframe
#' @returns
#' data.frame with:
#' mean of values per values_factors and/or weight_factors,
#' variance of values per values_factors and/or weight_factors,
#' skewness of values per values_factors and/or weight_factors,
#' kurtosis of values per values_factors and/or weight_factors,
#' @export
#' @import data.table
#' @examples
#' head(values)
#' Moments(
#' values = values$Values,
#' individual_id = NULL,
#' weight = NULL,
#' values_factors = NULL,
#' weight_factors = NULL
#' )

Moments <- function(
    values,
    individual_id = NULL,
    weight = NULL,
    values_factors = NULL,
    weight_factors = NULL
){
  # test data input validity
  if (is.null(values)) stop("We need the values to compute the moments of the distributions!")
  if (is.null(weight) && !is.null(weight_factors)) stop("You must specify the weighting dataframe to weight the values!")
  if (!is.null(weight) && (is.null(weight_factors) || is.null(individual_id))) stop("You must specify the weighting factors and/or the ID of individuals to weight the values!")
  if (!is.null(individual_id) && (nrow(data.frame(individual_id)) != nrow(data.frame(values)))) stop("values and individual id must have the same length!")
  if (!is.null(values_factors) && (nrow(data.frame(values_factors)) != nrow(data.frame(values)))) stop("values and values factors must have the same length!")
  if (!is.null(weight) && any(is.na(weight)) == TRUE) stop("weight data frame should not have NA values!")
  if (!is.null(weight) && any(duplicated(weight_factors)) == TRUE) stop("weight_factors must be a unique combination of factors!")
  if ((!is.null(weight) && !is.data.frame(weight)) || (!is.null(weight_factors) && !is.data.frame(weight_factors)) || (!is.null(values_factors) && !is.data.frame(values_factors))) stop("'weight', 'values_factors' & 'weight_factors' must be a dataframe!")
  
  # intern function to compute moments
  compute_moments <- function(x, w = NULL) {
    if (!is.null(w)) {
      w <- w / sum(w)
      m <- sum(x * w)
      v <- sum((x - m)^2 * w)
      s <- sum(((x - m)^3 / v^(3/2)) * w)
      k <- sum(((x - m)^4 / v^2) * w)
    } else {
      m <- mean(x)
      v <- var(x)
      s <- mean((x - m)^3) / (v^(3/2))
      k <- mean((x - m)^4) / (v^2)
    }
    data.frame(mean = m, variance = v, skewness = s, kurtosis = k)
  }
  
  # Case 1: Compute moment for a list of values ("values") ----
  if (is.null(weight) && is.null(values_factors)){
    return(compute_moments(x = values[!is.na(values)]))
  }
  
  # Case 2: Compute moments of a list of values ("values") observed under different conditions ("values_factors") ----
  if (!is.null(values_factors) && is.null(weight)) {
    dt_values <-  data.table::data.table(values_factors, values)
    groups <- colnames(values_factors)
    dt_values <- dt_values[!is.na(values)]
    return(dt_values[, as.list(compute_moments(values)), by = groups])
  }
  
  # Case 3: Compute weighting ("weight" & "weight_factors") moments with or without variability in the observation of values by individuals ("values_factors") ---- 
  # Factor for computation
  if (!is.null(weight)) {
    groups <- colnames(weight_factors)
    weight[] <- lapply(weight, function(x) as.numeric(as.character(x)))
    dt_weight <- data.table::data.table(weight_factors, weight)
    dt_weight <- data.table::melt(
      dt_weight,
      id.vars = colnames(weight_factors),
      variable.name = "individual_id",
      value.name = "weight"
    )
    if (!is.null(values_factors)) {
      dt_values <- data.table::data.table(individual_id, values_factors, values)
      dt_compil <- dt_weight[dt_values, on = c("individual_id", groups), nomatch = 0, allow.cartesian = FALSE]
      dt_compil <- dt_compil[!is.na(values)]
      return(dt_compil[, as.list(compute_moments(x = values, w = weight)), by = groups])
    } else {
      dt_values <- data.table::data.table(individual_id, values)
      dt_compil <- dt_weight[dt_values, on = "individual_id", nomatch = 0, allow.cartesian = FALSE]
      dt_compil <- dt_compil[!is.na(values)]
      return(dt_compil[, as.list(compute_moments(x = values, w = weight)), by = groups])
    }
  }
}