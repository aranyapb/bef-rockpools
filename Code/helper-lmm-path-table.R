#' Extract and Rename Predictor Terms from a Linear Model
#'
#' This function filters the coefficient table of a model to include only specified predictors,
#' renames the term column to 'path', and appends the response variable name to each path.
#'
#' @param predictors A character vector of predictor variable names to retain.
#' @param response A character string naming the response variable to append to each path.
#' @param model A fitted linear model object (e.g., output from `lm()`).
#'
#' @return A tibble with filtered and modified coefficient terms.
#' @export
#'
#' @examples
#' predictors <- c("loggamma", "logdepth", "pc1", "pc2")
#' response <- "logbiomass"
#' model <- lm(logbiomass ~ loggamma + logdepth + pc1 + pc2, data = your_data)
#' extract_predictor_paths(predictors, response, model)
extract_predictor_paths <- function(predictors, response, model) {
  library(dplyr)
  broom.mixed::tidy(model) |>
    dplyr::filter(term %in% predictors) |>
    dplyr::rename(path = term) |>
    dplyr::mutate(path = paste0(path, "-", response)) |>
    dplyr::select(-group, -effect) |>
    dplyr::mutate(dplyr::across(.cols = c("estimate", "std.error", "statistic"), ~round(.x, 6)))
}
