specificity <- function(confusion_matrix) {
  return( confusion_matrix[2, 2] / 
            (confusion_matrix[2,2] + confusion_matrix[1, 2]))
}

accuracy <- function(confusion_matrix) {
  return( (confusion_matrix[1, 1] + confusion_matrix[2,2]) / 
            sum(confusion_matrix))
}

precision <- function(confusion_matrix) {
  return( confusion_matrix[1, 1] / 
            (confusion_matrix[1,1] + confusion_matrix[1, 2]))
}

sensitivity <- function(confusion_matrix) {
  return( confusion_matrix[1,1] / 
            (confusion_matrix[1,1] + confusion_matrix[2, 1]))
}

negative_predictive_value <- function(confusion_matrix) {
  return( confusion_matrix[2,2] / 
            (confusion_matrix[2,2] + confusion_matrix[2, 1]))
}
