#' Calculate correlation between observed and predicted values
#'
#' @param obs Observed values.
#' @param pred Predicted values.
#' @param method Correlation method, e.g. "pearson", "spearman", or "kendall"
#'
#' @return The correlation coefficient.
#' @export
correlation <- function(obs, pred, method = "pearson"){
  return( cor(pred, obs, method = method))
}
