#'
#' @title ci_lmm
#' 
#' @description function to calculate the 95% confidence interval for a 
#' linear mixed model with a random intercept. The function uses the formulas
#' for standard linear model confidence intervals to calculate the confidence
#' interval using the average intercept from the lmm.
#' 
#' references for the confidence interval formula:
#' https://stats.stackexchange.com/questions/585660/what-is-the-formula-for-prediction-interval-in-multivariate-case#:~:text=It%20is%20given%20by%20%CB%86,is%20the%20number%20of%20regressors
#' https://real-statistics.com/multiple-regression/confidence-and-prediction-intervals/
#' https://book.stat420.org/multiple-linear-regression.html
#' 
#' calculate the mean squared residuals
#' https://book.stat420.org/multiple-linear-regression.html
#' 
#' @param obs_data - dataset used to fit the model
#' @param pred_data - dataset with new values of the predictor variables
#' @param obs_resp - name of the observed respond variable
#' @param pred - name of the predicted response variable
#' 

ci_lmm <- function(obs_data, pred_data, obs_resp, pred) {
  
  # calculate the number of predictor variables and number of datapoints
  p <- ncol(pred_data)
  n <- nrow(obs_data)
  
  MSres <- sum( (obs_data[[obs_resp]] - obs_data[[pred]])^2 )/(n - p)
  
  # get the upper tvalue
  tupp <- qt(p = 0.975, df = (n - p))
  tlow <- qt(p = 0.025, df = (n - p))
  
  # get the design matrix: X
  X_obs <- as.matrix(obs_data[, names(pred_data)])
  X_obs <- cbind(1, X_obs)
  X_pred <- t(cbind(1, as.matrix(pred_data)))
  
  # create upper ci values
  CI_upp <- vector(length = ncol(X_pred))
  CI_low <- vector(length = ncol(X_pred))
  
  for(i in 1:ncol(X_pred)) {
    
    Xt <- matrix(data = X_pred[, i], nrow = 1, ncol = ncol(X_obs))
    X <- matrix(data = X_pred[, i], nrow = ncol(X_obs), ncol = 1)
    
    # get the square root term
    sqrt_term <- sqrt(Xt%*%solve((t(X_obs)%*%X_obs))%*%X)
    
    # convert to a vector
    sqrt_term <- as.vector(sqrt_term)
    
    CI_upp[i] <- (tupp)*MSres*sqrt_term
    CI_low[i] <- (tlow)*MSres*sqrt_term
    
  }
  
  # pull CIs into a data.frame
  CI_df <- dplyr::tibble(CI_low = CI_low,
                         CI_upp = CI_upp)
  
  return(CI_df)
  
}
