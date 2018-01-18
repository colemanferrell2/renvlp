testcoef.penv <- function(m, L, R, A) {
# User needs to supply L, R and A as matrices.  
  if (is.null(m$Gamma)) stop("beta is a zero matrix, no test is interesting.")
  
  m$beta1 <- as.matrix(m$beta1)
  a <- dim(m$beta1)
  r <- a[1]
  p1 <- a[2]

  if (ncol(L) != r) stop("The size of L is not supported")
  if (nrow(R) != p1) stop("The size of R is not supported")
  if (nrow(L) != nrow(A) | ncol(R) != ncol(A)) stop("The size of A is not supported")
  
  tmp1 <- kronecker(t(R), L)
  Sigma <- tmp1 %*% tcrossprod(m$covMatrix[1:(r * p1), 1:(r * p1)], tmp1) / m$n
  tmp2 <- matrix(c(L %*% m$beta1 %*% R - A), nrow = 1)
  
  chisqStatistic <- tmp2 %*% tcrossprod(chol2inv(chol(Sigma)), tmp2)
  dof <- nrow(L) * ncol(R)
  pValue <- stats::pchisq(chisqStatistic, dof, lower.tail = F)
  covMatrix <- Sigma
  
  return(list(chisqStatistic = chisqStatistic, dof = dof, pValue = pValue, covMatrix = covMatrix))
}



