cv.genv <- function(X, Y, Z, u, m, nperm){

  groupind <- unique(Z)
  ZZ <- as.factor(Z)
  Z <- match(ZZ, levels(ZZ))
  ZZ <- as.factor(Z)
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(X)
  p <- nlevels(ZZ)
  
  if(u < 0 | u > r) {
    print("u should be an interger between 0 and r")
    skip <- 1
  } else {
    skip <- 0
  }
  if (n <= c) {
    print("This works when p < n")
    skip <- 1
  } else {
    skip <- 0
  }
  
  if(skip == 0) {
    
    prederr <- matrix(rep(0, m * nperm), ncol = nperm)
    PE <- rep(0, nperm)
    for (i in 1 : nperm) {
      id <- sample(n, n)
      Zn <- Z[id]
      Yn <- Y[id, ]
      Xn <- X[id, ]
      for (j in 1 : m) {
        id.test <- (floor((j - 1) * n / m) + 1) : ceiling(j * n / m)
        id.train <- setdiff(1:n, id.test)
        Z.train <- Zn[id.train]
        Y.train <- Yn[id.train, ]
        X.train <- Xn[id.train, ]
        Z.test <- Zn[id.test]
        Y.test <- Yn[id.test, ]
        X.test <- Xn[id.test, ]
        n.test <- length(id.test)
        
        fit <- genv(X.train, Y.train, Z.train, u, asy = F, fit = F)
        betahat <- fit$beta
        muhat <- fit$mu
        
        testn <- length(id.test)
        traceres <- 0
        for(l in 1 : testn){
          Ind <- fit$groupInd
          iG <- which(Ind == intersect(Z[l], Ind))
          resi <- Y.test[l, ] - t(muhat[ ,iG]) - X.test[l, ] %*% t(betahat[[iG]])
          traceres <- traceres + resi %*% t(resi)
        }
        
        prederr[j, i] <- traceres
      }
      PE[i] <- sqrt(sum(prederr[ ,i])/n)
    }
    
    out <- mean(PE)
    return(out)
  }
  
}