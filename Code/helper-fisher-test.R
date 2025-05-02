#' Fisher's Combined Probability Test
#'
#' Combines multiple independent p-values into a single p-value using Fisher's method.
#'
#' @param p_values A numeric vector of p-values (each between 0 and 1).
#'
#' @return A list containing:
#' \describe{
#'   \item{statistic}{The Fisher test statistic (-2 * sum(log(p)))} 
#'   \item{df}{Degrees of freedom (2 * number of p-values)}
#'   \item{p_value}{Combined p-value from the chi-squared distribution}
#' }
#' @export
#'
#' @examples
#' fisher_combined_test(c(0.01, 0.05, 0.10))
fisher_combined_test <- function(p_values) {
  if (any(p_values < 0 | p_values > 1)) {
    stop("All p-values must be in the interval (0, 1].")
  }
  
  # add tiny offset if p-value estimated at 0 due to numerical limitations
  p_value_trans <- ifelse(p_values == 0, p_values + 0.000001, p_values)
  
  fisher_stat <- -2 * sum(log(p_value_trans))
  df <- 2 * length(p_value_trans)
  p_value <- pchisq(fisher_stat, df = df, lower.tail = FALSE)
  
  return(list(
    statistic = round(fisher_stat, 4),
    df = df,
    p_value = round(p_value, 5)
  ))
}
