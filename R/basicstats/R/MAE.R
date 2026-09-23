#' Calculate mean absolute error
#'
#' @param obs Observed values.
#' @param pred Predicted values.
#'
#' @return The mean absolute error.
#' @export
MAE <- function(obs, pred){
  return(mean(abs(obs-pred)))
}

