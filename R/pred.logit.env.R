pred.logit.env <- function(m, Xnew) {
  if (is.null(m$ratio)) 
    stop("The asymptotic distribution part of logit.env is missing. Rerun logit.env with asy = T")
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
    value <- exp(t) / (1 + exp(t))
    covMatrix.estm <- NULL
  }
  else {
    Xnew <- as.matrix(Xnew)
    if (nrow(Xnew) == 1) 
      Xnew <- t(Xnew)
    t <- m$mu + crossprod(m$beta, Xnew)
    value <- exp(t) / (1 + exp(t))
    temp <- kronecker(diag(r), Xnew)
    temp1 <- crossprod(temp, m$covMatrix) %*% temp/n
    temp2 <- exp(t) / (1 + exp(t))^2
    covMatrix.estm <- crossprod(temp2, temp1)  %*% temp2
  }
  SE.estm <- sqrt(diag(covMatrix.estm))
  return(list(value = value, covMatrix.estm = covMatrix.estm, 
              SE.estm = SE.estm))
}
