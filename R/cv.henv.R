cv.henv <- function(X, Y, u, m, nperm) {

  groupind <- unique(X)
  XX <- as.factor(X)
  X <- match(XX, levels(XX))
  XX <- as.factor(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  p <- nlevels(XX)
  if(u < 0 | u > r){
    stop("u should be an interger between 0 and r")
  } 
  prederr <- matrix(rep(0, m * nperm), ncol = nperm)
  PE <- rep(0, nperm)
  for (i in 1 : nperm) {
    id <- sample(n, n)
    Xn <- X[id]
    Yn <- Y[id, ]
    for (j in 1:m) {
      id.test <- (floor((j - 1) * n/m) + 1):ceiling(j * n/m)
      id.train <- setdiff(1:n, id.test)
      X.train <- Xn[id.train]
      Y.train <- Yn[id.train, ]
      X.test <- Xn[id.test]
      Y.test <- Yn[id.test, ]
      n.test <- length(id.test)
      fit <- henv(X.train, Y.train, u, asy = F, fit = F)
      betahat <- fit$beta
      muhat <- fit$mu
      testn <- length(id.test)
      traceres <- 0
      for (l in 1:testn) {
        Ind <- fit$groupInd
        t <- which(Ind == intersect(X[l], Ind))
        pred <- fit$mug[ , t]
        resi <- Y.test[l, ] - t(pred)
        traceres <- traceres + resi %*% t(resi)
      }
      prederr[i] <- prederr[i] + traceres
    }
  }
  return(sqrt(mean(prederr/n)))
}