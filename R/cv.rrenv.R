cv.rrenv <- function(X, Y, u, d, m, nperm) {
  
  X <- as.matrix(X)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  
  prederr <- rep(0, nperm)
  
  for (i in 1:nperm)	{
    id <- sample(n, n)
    Xn <- as.matrix(X[id, ])
    Yn <- Y[id, ]
    for (j in 1:m) {
      id.test <- (floor((j - 1) * n / m) + 1) : ceiling(j * n / m)
      id.train <- setdiff(1:n, id.test)
      X.train <- Xn[id.train, ]
      Y.train <- Yn[id.train, ]
      X.test <- Xn[id.test, ]
      Y.test <- Yn[id.test, ]
      n.test <- length(id.test)
      fit <- rrenv(X.train, Y.train, u, d, asy = F)
      betahat <- fit$beta
      muhat <- fit$mu
      resi <- as.matrix(Y.test - matrix(1, n.test, 1) %*% t(muhat) - as.matrix(X.test) %*% t(betahat))
      prederr[i] <- prederr[i] + sum(resi^2)
    }
  }
  
  return(sqrt(mean(prederr / n)))
  
}
