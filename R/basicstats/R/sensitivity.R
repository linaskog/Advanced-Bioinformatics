#' Calculate sensitivity
#'
#' @param confusion_matrix Confusion Matrix 
#'   \preformatted{
#'   TP  FP
#'   FN  TN
#'   }
#' @return The sensitivity.
#' @export
sensitivity <- function(confusion_matrix) {
  return( confusion_matrix[1,1] / 
            (confusion_matrix[1,1] + confusion_matrix[2, 1]))
}
