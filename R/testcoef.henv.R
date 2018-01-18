testcoef.henv <- function(m, L, R, A) {
  
  if (is.null(m$Gamma)) 
    stop("beta is a zero matrix, no test is interesting.")
  r <- ncol(m$Sigma[[1]])
  n <- m$n
  ng <- m$ng
  p <- length(ng)
  if (ncol(L) != r) 
    stop("The size of L is not supported")
  if (nrow(R) != p) 
    stop("The size of R is not supported")
  if (nrow(L) != nrow(A)) 
    stop("The size of A is not supported")
  mm <- as.vector(kronecker(1 / sqrt(ng), matrix(1, r, 1)))
  multiplier <- diag(mm)
  tmp1 <- kronecker(t(R), L)
  aa <- r + 1
  bb <- r * (p + 1)
  cov <- m$covMatrix[aa : bb, aa : bb]
  Sigma <- tmp1 %*% multiplier %*% cov %*% multiplier %*% t(tmp1)
  
  tmp2 <- matrix(c(L %*% m$beta %*% R - A), nrow = 1)
  chisqStatistic <- tmp2 %*% tcrossprod(chol2inv(chol(Sigma)), tmp2)
  dof <- nrow(L) * ncol(R)
  pValue <- stats::pchisq(chisqStatistic, dof, lower.tail = F)
  covMatrix <- Sigma
  return(list(chisqStatistic = chisqStatistic, dof = dof, pValue = pValue, 
              covMatrix = covMatrix))
}