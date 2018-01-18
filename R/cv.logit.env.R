cv.logit.env <- function(X, Y, u, m, nperm){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- ncol(X)
  prederr <- rep(0, nperm)
  for (i in 1:nperm) {
    id <- sample(n, n)
    Xn <- X[id, ]
    Yn <- Y[id, ]
    Yn <- as.matrix(Yn)
    for (j in 1:m) {
      id.test <- (floor((j - 1) * n/m) + 1):ceiling(j * n/m)
      id.train <- setdiff(1:n, id.test)
      X.train <- Xn[id.train, ]
      Y.train <- Yn[id.train, ]
      X.test <- Xn[id.test, ]
      Y.test <- Yn[id.test, ]
      n.test <- length(id.test)
      fit <- logit.env(X.train, Y.train, u, asy = F)
      betahat <- fit$beta
      muhat <- fit$mu
      t.test <- matrix(1, n.test, 1) %*% t(muhat) + as.matrix(X.test) %*% betahat
      fit.test <- exp(t.test) / (1 + exp(t.test))
      resi <- as.matrix(Y.test - fit.test)
      sprederr <- apply(resi, 1, function(x) sum(x^2))
      prederr[i] <- prederr[i] + sum(sprederr)
    }
  }
  return(sqrt(mean(prederr/n)))
}