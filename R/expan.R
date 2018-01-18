expan <- function(d) {
  E <- matrix(0, d ^ 2, d * (d + 1) / 2)
  for (i in 1:d) {
    for (j in 1:d) {
      if (j >= i) {
        E[d * (i - 1) + j, (2 * d - i) / 2 * (i - 1) + j] <- 1
      } else {
        E[d * (i - 1) + j, (2 * d - j) / 2 * (j - 1) + i] <- 1
      }
    }
  }
  return(E)
}

