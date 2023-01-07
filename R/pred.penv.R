pred.penv <- function(m, X1new, X2new) {
  
  if(is.null(m$ratio)) stop("The asymptotic distribution part of penv is missing. Rerun penv with asy = T")
  
  r <- ncol(m$Sigma)
  n <- m$n
  if (is.null(m$Gamma)) {
    u <- 0
  } else { 
    u <- ncol(m$Gamma)
  }
  
  if (u == 0) {
    value <- m$mu + m$beta2 %*% X2new
    temp <- kronecker(X2new, diag(r)) 
    temp2 <- crossprod(temp, m$covMatrix) %*% temp / n
    covMatrix.estm <- m$Sigma / n + temp2
    covMatrix.pred <- (1 + 1 / n) * m$Sigma + temp2
  } else {
    X1new <- as.matrix(X1new)
    X2new <- as.matrix(X2new)
    if (nrow(X1new) == 1) X1new <- t(X1new)
    if (nrow(X2new) == 1) X2new <- t(X2new)
    Xnew <- rbind(X1new, X2new)
    value <- m$mu + m$beta1 %*% X1new + m$beta2 %*% X2new
    temp <- kronecker(Xnew, diag(r)) 
    covMatrix.estm <- m$Sigma / n + crossprod(temp, m$covMatrix) %*% temp / n
    covMatrix.pred <- (1 + 1 / n) * m$Sigma + crossprod(temp, m$covMatrix) %*% temp / n
  }
  
  SE.estm <- sqrt(diag(covMatrix.estm))
  SE.pred <- sqrt(diag(covMatrix.pred))
  return(
    list(
      value = value, covMatrix.estm = covMatrix.estm, SE.estm = SE.estm, covMatrix.pred = covMatrix.pred, SE.pred = SE.pred
    )
  )  
}
