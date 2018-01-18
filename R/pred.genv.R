pred.genv <- function(m, Xnew, Znew){
  if (is.null(m$ratio)) 
    stop("The asymptotic distribution part of env is missing. Rerun env with asy = T")
  r <- ncol(m$Sigma[[1]])
  n <- m$n
  ng <- m$ng
  Ind <- m$groupInd
  if (is.numeric(Ind)) {
    Ind <- sort(Ind, index.return = T)$x
  }
  t <- which(Ind == intersect(Znew, Ind))
  
  if (is.null(m$Gamma)) {
    u <- 0
  }
  else {
    u <- ncol(m$Gamma)
  }
  if (u == 0) {
    value <- m$mu
    covMatrix.estm <- m$Sigma/n
    covMatrix.pred <- (1 + 1/n) * m$Sigma
  }
  else {
    Xnew <- as.matrix(Xnew)
    if (nrow(Xnew) == 1) 
      Xnew <- t(Xnew)
    mu <- m$mu
    beta <- m$beta
    value <- mu[ ,t] + beta[[t]] %*% Xnew
    temp <- kronecker(Xnew, diag(r))
    covMatrix.estm <- m$Sigma[[t]]/ng[[t]] + crossprod(temp, 
                                            m$covMatrix[[t]]) %*% temp/ng[[t]]
    covMatrix.pred <- (1 + 1/ng[[t]]) * m$Sigma[[t]] + crossprod(temp, 
                                            m$covMatrix[[t]]) %*% temp/ng[[t]]
  }
  SE.estm <- sqrt(diag(covMatrix.estm))
  SE.pred <- sqrt(diag(covMatrix.pred))
  return(list(value = value, covMatrix.estm = covMatrix.estm, 
              SE.estm = SE.estm, covMatrix.pred = covMatrix.pred, SE.pred = SE.pred))
  
}