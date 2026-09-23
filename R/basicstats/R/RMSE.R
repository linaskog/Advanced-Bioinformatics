#' Calculate root mean squared error
#'
#' @param obs Observed values.
#' @param pred Predicted values.
#'
#' @return The root mean squared error.
#' @export
RMSE <- function(obs, pred){
  return( sqrt(mean((obs-pred)^2)))
}
