#' Predict from an lmer model with intervals using merTools::predictInterval
#'
#' This function generates predicted values and 95% confidence intervals from an `lmer` model
#' using parametric simulation via `merTools::predictInterval`. It supports both marginal and
#' group-specific (conditional) predictions.
#'
#' @param model A fitted `lmer` model object.
#' @param data The original data used to fit the model.
#' @param predictors Character vector of predictor names used in the model.
#' @param vary Character vector of predictors to vary across their range (for plotting).
#' @param n_points Integer, number of values to generate across each varying predictor.
#' @param nsim Number of simulations for prediction intervals.
#' @param include_re Logical, whether to include random effects (conditional predictions).
#' @param group_levels Optional named list specifying group levels for which to generate predictions.
#' @param seed Random seed for reproducibility.
#'
#' @return A `data.frame` with predictor values and prediction intervals (fit, lwr, upr).
#' @export
predict_lmer_with_intervals_merTools <- function(model, data, predictors,
                                                 vary, n_points = 100, nsim = 1000,
                                                 include_re = TRUE,
                                                 group_levels = NULL,
                                                 seed = 123) {
  set.seed(seed)
  
  group_vars <- names(ranef(model))
  
  if (include_re && is.null(group_levels)) {
    stop("If include_re = TRUE, you must specify group_levels (e.g., list(Group = c('A', 'B')))")
  }
  
  if (include_re) {
    newdata_list <- list()
    
    for (i in seq_along(group_levels[[1]])) {
      newdata_tmp <- data.frame(matrix(ncol = length(predictors), nrow = n_points))
      names(newdata_tmp) <- predictors
      
      for (var in predictors) {
        if (var %in% vary) {
          newdata_tmp[[var]] <- seq(min(data[[var]], na.rm = TRUE),
                                    max(data[[var]], na.rm = TRUE),
                                    length.out = n_points)
        } else {
          newdata_tmp[[var]] <- mean(data[[var]], na.rm = TRUE)
        }
      }
      
      for (g in group_vars) {
        newdata_tmp[[g]] <- group_levels[[g]][i]
      }
      
      newdata_tmp$.group_id <- group_levels[[g]][i]
      newdata_list[[i]] <- newdata_tmp
    }
    
    newdata <- do.call(rbind, newdata_list)
    
    # Ensure group vars are factors with correct levels
    for (h in group_vars) {
      newdata[[h]] <- factor(newdata[[h]], levels = levels(data[[h]]))
    }
    
  } else {
    newdata <- data.frame(matrix(ncol = length(predictors), nrow = n_points))
    names(newdata) <- predictors
    
    for (var in predictors) {
      if (var %in% vary) {
        newdata[[var]] <- seq(min(data[[var]], na.rm = TRUE),
                              max(data[[var]], na.rm = TRUE),
                              length.out = n_points)
      } else {
        newdata[[var]] <- mean(data[[var]], na.rm = TRUE)
      }
    }
    
    for (g in group_vars) {
      newdata[[g]] <- sample(x = data[[g]], 1)
    }
  }
  
  # Run merTools::predictInterval
  pred_mod <- merTools::predictInterval(
    merMod = model,
    newdata = newdata,
    level = 0.95,
    n.sims = nsim,
    stat = "mean",
    include.resid.var = FALSE,
    which = if (include_re) "full" else "fixed"
  )
  
  newdata$predicted_mean <- pred_mod$fit
  newdata$lower <- pred_mod$lwr
  newdata$upper <- pred_mod$upr
  
  return(newdata)
}
