pred.rrenv <- function(m, Xnew) {
  if (is.null(m$ratio))
    stop("The asymptotic distribution part of rrenv is missing. Rerun rrenv with asy = T")
  
  r <- ncol(m$Sigma)
  n <- m$n
  if (is.null(m$Gamma)) {
    u <- 0
  } else {
    u <- ncol(m$Gamma)
  }
  
  if (u == 0) {
    value <- m$mu
    covMatrix.estm <- m$Sigma / n
    covMatrix.pred <- (1 + 1 / n) * m$Sigma
    
  } else {
    Xnew <- as.matrix(Xnew)
    if (nrow(Xnew) == 1)
      Xnew <- t(Xnew)
    value <- m$mu + m$beta %*% Xnew
    temp <- kronecker(Xnew, diag(r))
    covMatrix.estm <-
      m$Sigma / n + crossprod(temp, m$covMatrix) %*% temp / n
    covMatrix.pred <-
      (1 + 1 / n) * m$Sigma + crossprod(temp, m$covMatrix) %*% temp / n
  }
  
  SE.estm <- sqrt(diag(covMatrix.estm))
  SE.pred <- sqrt(diag(covMatrix.pred))
  return(
    list(
      value = value, covMatrix.estm = covMatrix.estm, SE.estm = SE.estm, covMatrix.pred = covMatrix.pred, SE.pred = SE.pred
    )
  )
  
}
