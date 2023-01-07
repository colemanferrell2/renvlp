pred.pois.env <- function(m, Xnew) {
  if (is.null(m$ratio)) 
    stop("The asymptotic distribution part of pois.env is missing. Rerun pois.env with asy = T")
  n <- m$n
  r <- 1
  if (is.null(m$Gamma)) {
    u <- 0
  }
  else {
    u <- ncol(m$Gamma)
  }
  if (u == 0) {
    t <- m$mu
    value <- exp(t)
    covMatrix.estm <- NULL
  }
  else {
    Xnew <- as.matrix(Xnew)
    if (nrow(Xnew) == 1) 
      Xnew <- t(Xnew)
    t <- m$mu + crossprod(m$beta, Xnew)
    value <- exp(t)
    temp <- kronecker(diag(r), Xnew)
    temp1 <- crossprod(temp, m$covMatrix) %*% temp/n
    covMatrix.estm <- crossprod(exp(t), temp1)  %*% exp(t)
  }
  SE.estm <- sqrt(diag(covMatrix.estm))
  return(list(value = value, covMatrix.estm = covMatrix.estm, 
              SE.estm = SE.estm))
}
