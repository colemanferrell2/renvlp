testcoef.genv <- function(m, L, R, A){
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
  if (nrow(L) != nrow(A) | ncol(R) != ncol(A)) 
    stop("The size of A is not supported")
  tmp1 <- kronecker(t(R), L)
  
  Sigma = tmp2 = pValue = covMatrix = chisqStatistic <- list(length = p)
  for (i in 1 : p) {
    Sigma[[i]] <- tmp1 %*% tcrossprod(m$covMatrix[[i]], tmp1)/m$ng[[i]]
    tmp2[[i]] <- matrix(c(L %*% m$beta[[i]] %*% R - A), nrow = 1)
    chisqStatistic[[i]] <- tmp2[[i]] %*% tcrossprod(chol2inv(chol(Sigma[[i]])), 
                                        tmp2[[i]])
    dof <- nrow(L) * ncol(R)
    pValue[[i]] <- stats::pchisq(chisqStatistic[[i]], dof, lower.tail = F)
    covMatrix[[i]] <- Sigma[[i]]
  }
  return(list(chisqStatistic = chisqStatistic, dof = dof, pValue = pValue, 
              covMatrix = covMatrix))
}