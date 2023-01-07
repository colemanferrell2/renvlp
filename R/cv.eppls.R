cv.eppls <- function(X1, X2, Y, u, m, nperm) {
  
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p1 <- ncol(X1)
  p2 <- ncol(X2)
  
  prederr <- rep(0, nperm)

  for (i in 1:nperm)	{
    id <- sample(n, n)
    X1n <- as.matrix(X1[id, ])
    X2n <- as.matrix(X2[id, ])
    Yn <- Y[id, ]
    for (j in 1:m) {
      id.test <- (floor((j - 1) * n / m) + 1) : ceiling(j * n / m)
      id.train <- setdiff(1:n, id.test)
      X1.train <- X1n[id.train, ]
      X2.train <- X2n[id.train, ]
      Y.train <- Yn[id.train, ]
      X1.test <- X1n[id.test, ]
      X2.test <- X2n[id.test, ]
      Y.test <- Yn[id.test, ]
      n.test <- length(id.test)
      fit <- eppls(X1.train, X2.train, Y.train, u, asy = F)
      beta1hat <- fit$beta1
      beta2hat <- fit$beta2
      muhat <- fit$muY
      resi <- as.matrix(Y.test - matrix(1, n.test, 1) %*% t(muhat) - as.matrix(X1.test) %*% t(beta1hat) - as.matrix(X2.test) %*% t(beta2hat))
      sprederr <- apply(resi, 1, function(x) sum(x ^ 2))
      prederr[i] <- prederr[i] + sum(sprederr)
    }
  }
  
  return(sqrt(mean(prederr / n)))
  
}