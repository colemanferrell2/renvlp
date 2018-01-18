pred.henv <- function (m, Xnew) {
  
  if(is.null(m$ratio)) stop("The asymptotic distribution part of env is missing. Rerun henv with asy = T")
  
  r <- ncol(m$Sigma[[1]])
  n <- m$n
  ng <- m$ng
  Ind <- m$groupInd
  if (is.numeric(Ind)) {
      Ind <- sort(Ind, index.return = T)$x
  }
  t <- which(Ind == intersect(Xnew, Ind))
  covMatrix.estm <- NULL
  covMatrix.pred <- NULL
  SE.estm <- NULL
  SE.pred <- NULL
  if (is.null(m$Gamma)) {
    u <- 0
  }
  else {
    u <- ncol(m$Gamma)
  }
  if (u == 0) {
    value <- m$mu
    Sigma <- m$covMatrix[1 : r, 1 : r]
    covMatrix.estm <- Sigma / n
    covMatrix.pred <- Sigma / n + m$Sigma[[t]]
    SE.estm <- sqrt(diag(covMatrix.estm))
    SE.pred <- sqrt(diag(covMatrix.pred))
  } else {
    value <- m$mug[ , t]
    aa = t * r + 1
    bb = (t + 1) * r
    covMatrix.estm <- (m$covMatrix[1 : r, 1 : r] + m$covMatrix[aa : bb, aa : bb] 
                       + m$covMatrix[1 : r, aa : bb] 
                       + m$covMatrix[aa : bb, 1 : r] / ng[t])
    covMatrix.pred <- (m$covMatrix[1 : r, 1 : r] + m$covMatrix[aa : bb, aa : bb] 
                       + m$covMatrix[1 : r, aa : bb] 
                       + m$covMatrix[aa : bb, 1 : r] / ng[t] 
                       + m$Sigma[[t]])
    SE.estm <- sqrt(diag(covMatrix.estm))
    SE.pred <- sqrt(diag(covMatrix.pred))
  }  
  return(list(value = value, 
              covMatrix.estm = covMatrix.estm, SE.estm = SE.estm, 
              covMatrix.pred = covMatrix.pred, SE.pred = SE.pred))
}