
#' Fit and Compare Linear Mixed Models with DHARMa Diagnostics
#'
#' This function fits multiple linear mixed models to the same dataset using a list of model formulas.
#' It evaluates model assumptions using DHARMa, compares models by AIC, and returns a summary table.
#'
#' @param formulas A named list of model formulas.
#' @param data The dataset used for model fitting.
#' @param n_sim Number of simulations for DHARMa residual diagnostics. Default is 1000.
#'
#' @return A list with elements:
#' \describe{
#'   \item{models}{A named list of fitted lmer models}
#'   \item{residuals}{A list of DHARMa simulation objects}
#'   \item{comparison}{A data frame comparing AIC, singularity, and assumption test results}
#' }
#' @export
fit_and_compare_lmm_models <- function(formulas, data, n_sim = 1000) {
  stopifnot(requireNamespace("lme4"), requireNamespace("DHARMa"), requireNamespace("dplyr"))
  
  model_fits <- list()
  model_checks <- list()
  
  for (name in names(formulas)) {
    model_fits[[name]] <- lme4::lmer(
      formula = formulas[[name]],
      data = data,
      REML = FALSE
    )
    model_checks[[name]] <- DHARMa::simulateResiduals(fittedModel = model_fits[[name]], n = n_sim)
  }
  
  model_aics <- sapply(model_fits, AIC)
  
  model_comparison <- dplyr::tibble(
    model = names(model_aics),
    AIC = as.numeric(model_aics),
    singular = sapply(model_fits, lme4::isSingular)
  )
  
  # Run DHARMa tests
  dharma_tests <- lapply(model_checks, function(x) {
    quant_test <- DHARMa::testQuantiles(x, plot = FALSE)
    res_test <- DHARMa::testResiduals(x, plot = FALSE)
    
    dplyr::tibble(
      quantile_test = quant_test$p.value,
      uniformity_test = res_test$uniformity$p.value,
      dispersion_test = res_test$dispersion$p.value,
      outlier_test = res_test$outliers$p.value
    )
  })
  
  dharma_tests_df <- dplyr::bind_cols(
    dplyr::tibble(model = names(dharma_tests)),
    dplyr::bind_rows(dharma_tests)
  )
  
  model_comparison <- dplyr::full_join(model_comparison, dharma_tests_df, by = "model")
  
  # Logical flag for whether all assumption tests passed
  model_comparison$assumptions_met <- with(
    model_comparison,
    ifelse(
      quantile_test > 0.05 &
        uniformity_test > 0.05 &
        dispersion_test > 0.05 &
        outlier_test > 0.05,
      "yes", "no"
    )
  )
  
  # Order by AIC
  model_comparison <- model_comparison[order(model_comparison$AIC), ]
  
  return(list(
    models = model_fits,
    residuals = model_checks,
    comparison = model_comparison
  ))
}
