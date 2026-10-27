#' Calculate precision
#'
#' @param confusion_matrix Confusion Matrix 
#'   \preformatted{
#'   TP  FP
#'   FN  TN
#'   }
#' @return The precision.
#' @export
precision <- function(confusion_matrix) {
  return( confusion_matrix[1, 1] / 
            (confusion_matrix[1,1] + confusion_matrix[1, 2]))
}
