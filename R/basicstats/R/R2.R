#' Calculate R-squared
#'
#' @param obs Observed values.
#' @param pred Predicted values.
#'
#' @return The R-squared value.
#' @export
R2 <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}
