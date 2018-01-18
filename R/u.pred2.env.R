u.pred2.env <- function(X, Y, Xnew, alpha = 0.01) {
  
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

 u <- u.penv(X1, X2, Y)
  
 return(list(u.aic = u$u.aic, u.bic = u$u.bic, u.lrt = u$u.lrt, loglik.seq = u$loglik.seq, aic.seq = u$aic.seq, bic.seq = u$bic.seq))
}
