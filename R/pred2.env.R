pred2.env <- function(X, Y, u, Xnew) {
  
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
  
  # u <- u.penv(X1, X2, Y)$u.bic
  fit <- penv(X1, X2, Y, u)
  
  X1new <- Ainv[1, ] %*% Xnew
  X2new <- Ainv[2:p, ] %*% Xnew
  tmp <- pred.penv(fit, X1new, X2new)
  
  return(list(value = tmp$value, covMatrix.estm = tmp$covMatrix.estm, SE.estm = tmp$SE.estm, covMatrix.pred = tmp$covMatrix.pred, SE.pred = tmp$SE.pred))
  
}
