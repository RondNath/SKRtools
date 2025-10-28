#' @title SKR analysis - compute regression parameters of the SKR
#' @description Launch the SKR analysis of datasets with skewness and kurtosis values, based on explanatory factors (e.g., spatial, temporal, context, ...)
#' @param values vector of values observed/measured per individual (length format)
#' @param individual_id vector of Individual name (length format, same length as values, must be completed if weight =! NULL and correspond to the column name in the weighting data frame)
#' @param weight dataframe of the weight to attribute to each individual for each factors of observation (e.g., frequency, abundance, ..., NO NA) (wide format, column name = Individual observed, must correspond to the row name in the values data frame)
#' @param values_factors dataframe of the factors of observations of the values vector
#' @param weight_factors dataframe of the factors of observations of the weight dataframe
#' @param workers Cores used for parallelisation (used 1 for no parallelisation)
#' @returns
#' data.frame with:
#' mean of values per values_factors and/or weight_factors,
#' variance of values per values_factors and/or weight_factors,
#' skewness of values per values_factors and/or weight_factors,
#' kurtosis of values per values_factors and/or weight_factors,
#' @export
#' @importFrom dplyr filter select left_join %>%
#' @importFrom tidyr unite
#' @importFrom future plan
#' @importFrom parallel detectCores mclapply
#' @importFrom stats mean var na.omit
#' @examples
#' head(values)
#' Moments(
#' values = values$Values,
#' individual_id = NULL,
#' weight = NULL,
#' values_factors = NULL,
#' weight_factors = NULL,
#' workers = 1
#' )

Moments <- function(
    values,
    individual_id = NULL,
    weight = NULL,
    values_factors = NULL,
    weight_factors = NULL,
    workers = parallel::detectCores() - 1
){
  # test data input validity
  if (is.null(values)){
    stop("We need the values to compute the moments of the distributions!")
  }
  if (is.null(weight) && !is.null(weight_factors)){
    stop("You must specify the weighting dataframe to weight the values!")
  }
  if (!is.null(weight) && (is.null(weight_factors) || is.null(individual_id))){
    stop("You must specify the weighting factors and/or the ID of individuals to weight the values!")
  }
  if (!is.null(individual_id) && (nrow(data.frame(individual_id)) != nrow(data.frame(values)))){
    stop("values and individual id must have the same length!")
  }
  if (!is.null(values_factors) && (nrow(data.frame(values_factors)) != nrow(data.frame(values)))){
    stop("values and values factors must have the same length!")
  }
  if (!is.null(weight) && any(is.na(weight)) == TRUE){
    stop("weight data frame should not have NA values!")
  }
  if (!is.null(weight) && any(duplicated(weight_factors)) == TRUE){
    stop("weight_factors must be a unique combination of factors!")
  }
  
  # session type for parallelisation
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, 
                 workers = workers)
  } else {
    future::plan(future::multicore, 
                 workers = workers)
  }
  
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
  if (is.null(weight) && is.null(values_factors)) {
    # remove NA
    values <- na.omit(values) 
    # compute moments
    return(compute_moments(x = values))
  }
  
  # Case 2: Compute moments of a list of values ("values") observed under different conditions ("values_factors") ----
  if (!is.null(values_factors) && is.null(weight)) {
    # remove NA
    df_values <- cbind(data.frame(values_factors) %>% 
                         tidyr::unite("Factors", 
                                      dplyr::everything(), 
                                      sep = "_", 
                                      remove = FALSE), 
                       values) %>% 
      na.omit()
    
    # compute moments per factors
    group_factor <- unique(data.frame(values_factors) %>% 
                             tidyr::unite("Factors", 
                                          dplyr::everything(), 
                                          sep = "_", 
                                          remove = FALSE))
    results <- parallel::mclapply(group_factor$Factors, 
                        function(g) {
                          compute_moments(df_values$values[df_values$Factors == g])
                        }, 
                        mc.cores = workers)
    
    return(data.frame(group_factor %>% 
                        dplyr::select(-Factors), 
                      do.call(rbind, 
                              results)))
  }
  
  # Case 3: Compute weighting ("weight" & "weight_factors") moments with or without variability in the observation of values by individuals ("values_factors") ---- 
  # Factor for computation
  if (!is.null(weight)) {
    df_weight <- cbind(data.frame(weight_factors) %>%
                         tidyr::unite("Factors", 
                               dplyr::everything(), 
                               sep = "_", 
                               remove = FALSE), 
                       weight)
    if (is.null(values_factors)) {
      df_values <- data.frame(individual_id, 
                              values)
    } else {
      df_values <- data.frame(data.frame(values_factors) %>%
                                tidyr::unite("Factors", 
                                      dplyr::everything(), 
                                      sep = "_", 
                                      remove = FALSE), 
                              individual_id, 
                              values)
    }
    # check if values and weight are linked by the same factors
    if (!is.null(values_factors)) {
      if (setequal(df_weight$Factors, 
                   unique(df_values$Factors)) == FALSE) {
        stop("weight_factors and values_factors must be the same!")
      }
    }
    group_factor <- data.frame(weight_factors) %>% 
                             tidyr::unite("Factors", 
                                          dplyr::everything(), 
                                          sep = "_", 
                                          remove = FALSE)
    # compil weight and values per Factors
    results <- parallel::mclapply(group_factor$Factors, 
                        function(g) {
      w <- df_weight %>% 
        dplyr::filter(Factors == g)
      wv <- data.frame(weight = t(w %>% 
                                    dplyr::select(-colnames(group_factor))))
      wv$individual_id <- rownames(wv)
      
      if (!is.null(values_factors)) {
        v <- df_values %>%
          dplyr::filter(Factors == g) %>%
          dplyr::select(-colnames(group_factor))
      } else {
        v <- df_values
      }
      # remove NA
      merged_wv <- dplyr::left_join(v, 
                          wv, 
                          by = "individual_id") %>% 
        na.omit()
      # compute moments
      compute_moments(x = merged_wv$values, 
                      w = merged_wv$weight)
    }, 
    mc.cores = workers)
    
    return(data.frame(group_factor %>% 
                        dplyr::select(-Factors), 
                      do.call(rbind, 
                              results)))
  }
}
