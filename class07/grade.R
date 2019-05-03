# Author: Akshara Balachandra
# Date: 04/24/19

grade <- function(grades) {
  
  min.grade.ind <- ifelse(sum(is.na(grades)), # if there is an NA, find that NA
      which(is.na(grades)),
      which(grades == min(grades))[1]) # otherwise, find the minimum grade
  
  keep <- rep(T, length(grades))
  keep[min.grade.ind] <- F
  return(grades[keep])
}