boot.logit.env <- function(X, Y, u, B) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  bootenv <- function(i) {
    sam <- sample(1:n, n, replace = T)
    Y.boot <- Y[sam]
    X.boot <- X[sam, ]
    return(c(logit.env(X.boot, Y.boot, u, asy = F)$beta))
  }
  bootbeta <- lapply(1:B, function(i) bootenv(i))
  bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)
  bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = p)
  return(bootse)
}