contr <- function(d) {
  C <- matrix(0, d * (d + 1) / 2, d ^ 2)
  for (i in 1:d) {
    for (j in 1:d) {
      if (j == i) {
        C[(2 * d - i) / 2 * (i - 1) + j, d * (i - 1) + j] <- 1
      } else if (j > i) {
        C[(2 * d - i) / 2 * (i - 1) + j, d * (i - 1) + j] <- 1 / 2
      } else {
        C[(2 * d - j) / 2 * (j - 1) + i, d * (i - 1) + j] <- 1 / 2
      }
    }
  }
  return(C)  
}
