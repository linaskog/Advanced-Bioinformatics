#' Calculate Negative Predictive Value
#'
#' @param confusion_matrix Confusion Matrix 
#'   \preformatted{
#'   TP  FP
#'   FN  TN
#'   }
#' @return The Negative Predictive Value.
#' @export
NPV <- function(confusion_matrix) {
  return( confusion_matrix[2,2] / 
            (confusion_matrix[2,2] + confusion_matrix[2, 1]))
}