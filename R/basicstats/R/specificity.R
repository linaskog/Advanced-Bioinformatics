#' Calculate specificity
#'
#' @param confusion_matrix Confusion Matrix 
#'   \preformatted{
#'   TP  FP
#'   FN  TN
#'   }
#' @return The specificity.
#' @export
specificity <- function(confusion_matrix) {
  return( confusion_matrix[2, 2] / 
            (confusion_matrix[2,2] + confusion_matrix[1, 2]))
}


