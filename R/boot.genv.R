boot.genv <- function(X, Y, Z, u, B) {
  
  ZZ <- as.factor(Z)
  a <- dim(Y)
  n <- a[1]
  r <- a[2]
  c <- ncol(X)
  p <- nlevels(ZZ)
  
  if(u < 0 | u > r){
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
    
    ncumx <- c()
    for (i in 1 : p) {
      ncumx[i] <- length(which(ZZ == as.numeric(levels(ZZ)[i])))
    }
    ncum <- cumsum(ncumx)
    ng <- diff(c(0, ncum))
    sortz <- sort(Z, index.return = T)
    Zs <- sortz$x
    ind <- sortz$ix
    
    fit <- genv(X, Y, Z, u, asy = F, fit = T)
    Yfit <- fit$Yfit
    res <- Y - Yfit
    
    bootgenv <- function(i) {
      out <- list(length = p)
      res.boot <- matrix(rep(0, n * r), ncol = r)
      for (j in 1 : p) {
        if(j > 1) {
          res.boot[ind[(ncum[j - 1] + 1):ncum[j]], ] <- as.matrix(res[sample(ind[(ncum[j - 1] + 1):ncum[j]], ng[j], replace = T), ], ncol = r)
        } else {
          res.boot[ind[1:ncum[1]], ] <- as.matrix(res[sample(ind[1:ncum[1]], ng[1], replace = T), ], ncol = r)
        }
      }
      Y.boot <- Yfit + res.boot
      for(k in 1 : p) {
        out[[k]] <- genv(X, Y.boot, Z, u, asy = F, fit = T)$beta[[k]]
      }
      return(out)
    }
    
    bootsebeta <- list(length = p)
    for (k in 1 : p) {
      out1 <- lapply(1 : B, function(i) bootgenv(i)[[k]])
      bootbeta <- matrix(unlist(out1), nrow = B, byrow = TRUE)
      bootsebeta[[k]] <- matrix(apply(bootbeta, 2, stats::sd), nrow = r)
    }
    
    return(list(bootse = bootsebeta))
  }
  
}