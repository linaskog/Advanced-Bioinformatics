#' Calculate accuracy
#'
#' @param confusion_matrix Confusion Matrix 
#'   \preformatted{
#'   TP  FP
#'   FN  TN
#'   }
#' @return The accuracy.
#' @export
accuracy <- function(confusion_matrix) {
  return( (confusion_matrix[1, 1] + confusion_matrix[2,2]) / 
            sum(confusion_matrix))
}
