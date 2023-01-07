pred.sxenv <- function (m, Xnew) 
{
  if (is.null(m$ratio)) 
    stop("The asymptotic distribution part of sxenv is missing. Rerun sxenv with asy = T")
  n <- m$n
  r <- ncol(m$SigmaYcX)
  if (is.null(m$Gamma)) {
    u <- 0
  }
  else {
    u <- ncol(m$Gamma)
  }
  if (u == 0) {
    value <- m$muY
    covMatrix.estm <- m$SigmaYcX/n
    covMatrix.pred <- (1 + 1/n) * m$SigmaYcX
  }
  else {
    Xnew <- as.matrix(Xnew)
    if (nrow(Xnew) == 1) 
      Xnew <- t(Xnew)
    value <- m$muY - crossprod(m$beta, m$muX) + crossprod(m$beta, Xnew)
    temp <- kronecker(diag(r), Xnew)
    covMatrix.estm <- m$SigmaYcX/n + crossprod(temp, m$covMatrix) %*% 
      temp/n
    covMatrix.pred <- (1 + 1/n) * m$SigmaYcX + crossprod(temp, 
                                                         m$covMatrix) %*% temp/n
  }
  SE.estm <- sqrt(diag(covMatrix.estm))
  SE.pred <- sqrt(diag(covMatrix.pred))
  return(list(value = value, covMatrix.estm = covMatrix.estm, 
              SE.estm = SE.estm, covMatrix.pred = covMatrix.pred, SE.pred = SE.pred))
}
