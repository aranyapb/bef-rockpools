
#' Plot DFBETAs for lmer Model (Leave-One-Out Influence Diagnostics)
#'
#' This function performs a leave-one-out diagnostic for a fitted `lmer` model,
#' computing the change in fixed effect estimates (DFBETAs) caused by removing each observation.
#' It produces a plot for each fixed effect showing the standardized influence of each observation.
#'
#' @param model A fitted `lmer` model object from the `lme4` package.
#' @param data The original dataset used to fit the model.
#' @param formula An optional formula object. If not provided, the formula is extracted from `model`.
#' @param id_var Optional character. The name of a column in `data` to label or track observations. Currently not used for plotting but can be included for future features.
#'
#' @return Invisibly returns a matrix of DFBETAs with rows corresponding to observations and columns to fixed effects.
#'
#' @details
#' DFBETAs represent the standardized change in each fixed effect coefficient when a single observation is removed.
#' This is a useful diagnostic for detecting influential individual observations in linear mixed models.
#'
#' An observation is often considered influential if the absolute DFBETA exceeds 1.
#' The function will generate a plot for each fixed effect with horizontal reference lines at ±1.
#'
#' @importFrom lme4 lmer fixef
#' @export
#'
#' @examples
#' \dontrun{
#' library(lme4)
#' model <- lmer(logbiomass ~ logalpha + (1 | Inselberg), data = my_data)
#' dfb <- plot_lmer_dfbetas(model, data = my_data)
#' }
plot_lmer_dfbetas <- function(model, data, formula = NULL, id_var = NULL) {
  if (!requireNamespace("lme4", quietly = TRUE)) stop("lme4 package required.")
  
  # Extract model formula if not provided
  if (is.null(formula)) {
    formula <- formula(model)
  }
  
  n <- nrow(data)
  beta_mat <- matrix(NA, nrow = n, ncol = length(fixef(model)))
  colnames(beta_mat) <- names(fixef(model))
  
  message("Running leave-one-out loop for ", n, " observations...")
  
  for (i in 1:n) {
    data_i <- data[-i, ]
    model_i <- try(lme4::lmer(formula, data = data_i), silent = TRUE)
    
    if (inherits(model_i, "try-error")) {
      warning("Model failed to converge at observation ", i)
      next
    }
    
    beta_mat[i, ] <- fixef(model_i)
  }
  
  # Compute DFBETAs
  se <- summary(model)$coefficients[, "Std. Error"]
  dfbetas_mat <- scale(beta_mat, center = fixef(model), scale = se)
  
  # Plot DFBETAs
  op <- par(mfrow = c(1, ncol(dfbetas_mat)))
  on.exit(par(op))
  
  for (j in 1:ncol(dfbetas_mat)) {
    plot(dfbetas_mat[, j],
         main = paste("DFBETA for", colnames(dfbetas_mat)[j]),
         xlab = "Observation",
         ylab = "DFBETA",
         pch = 20)
    abline(h = c(-1, 1), lty = 2, col = "red")
  }
  
  invisible(dfbetas_mat)
}
