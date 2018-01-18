weighted.pred.env <- function(X, Y, Xnew) {
  
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  Xnew <- as.matrix(Xnew)
  if (nrow(Xnew) == 1) Xnew <- t(Xnew)
  A <- qr.Q(qr(Xnew), complete = TRUE)
  Ainv <- solve(A)
  Z <- tcrossprod(X, Ainv)
  X1 <- Z[, 1]
  X2 <- Z[, 2:p]
  
  fit <- weighted.penv(X1, X2, Y)
  X1new <- Ainv[1, ] %*% Xnew
  X2new <- Ainv[2:p, ] %*% Xnew
  value <- fit$mu + fit$beta1 %*% X1new + fit$beta2 %*% X2new

  return(value)
  
}
