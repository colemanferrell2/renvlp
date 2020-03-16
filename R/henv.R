henv <- function(X, Y, u, asy = TRUE, fit = TRUE, init = NULL) {
  
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
  
  if (!is.null(init)) {
      if (nrow(init) != r || ncol(init) != u) stop("The dimension of init is wrong.")
  }
  
  ncumx <- c()
  for (i in 1 : p) {
    ncumx[i] <- length(which(XX == as.numeric(levels(XX)[i])))
  }
  ncum <- cumsum(ncumx)
  ng <- diff(c(0, ncum))
  fracN <- ng / n
  sortx <- sort(X, index.return = T)
  Xs <- sortx$x
  ind <- sortx$ix
  Ys <- Y[ind, ]
  mY <- apply(Y, 2, mean)
  
  sigres = invres <- list(length = p)
  mYg <- matrix(rep(0, r*p), ncol = p)
  for (i in 1 : p) {
    if(i > 1) {
      yy <- Ys[(ncum[i-1] + 1):ncum[i], ]
      sigres[[i]] <- stats::cov(yy) * (ng[i] - 1)/ng[i]
      mYg[ , i] <- apply(yy, 2, mean)
      invres[[i]] <- chol2inv(chol(sigres[[i]]))
    } else {
      yy <- Ys[1:ncum[i], ]
      sigres[[i]]<- stats::cov(yy) * (ng[i] - 1)/ng[i]
      mYg[ , i] <- apply(yy, 2, mean)
      invres[[i]] <- chol2inv(chol(sigres[[i]]))
    }
  }
  sigY <- stats::cov(Y) * (n - 1)/n
  eigtemY <- eigen(sigY)$values
  logDetSigY <- log(prod(eigtemY[eigtemY > 0]))
  invsigY <- chol2inv(chol(sigY))
  
  U = M <- list(length = p)
  for (i in 1 : p){
    M[[i]] <- sigres[[i]]
    U[[i]] <- sigY - M[[i]]
  }
  MU <- sigY
  tmp <- henvMU(M, U, MU, u, n, ng, p, initial = init)
  
  
  Gammahat <- tmp$Gammahat
  Gamma0hat <- tmp$Gamma0hat
  covMatrix <- NULL
  asySE <- NULL
  ratio <- NULL
  Yfit <- NULL
  
  if (u == 0) {
    Sigmahat <- sigY
    etahat <- NULL
    Omegahat <- NULL
    Omega0hat <- sigY
    betahat <- matrix(0, r, p)
    muhat <- as.matrix(mY)
    mughat <- mY %*% matrix(1, 1, p)
    loglik <- -n * r * (1 + log(2 * pi)) / 2 - n * logDetSigY / 2
    if (asy == T) {
      ratio <- matrix(1, r, p)
    }
    if (fit == T) {
      Yfit <- matrix(0, n, r)
    }
  }
  else if (u == r) {
    Sigmahat <- sigres
    muhat <- as.matrix(mY)
    mughat <- mYg 
    betahat <- mughat - mY %*% matrix(1, 1, p)
    etahat <- betahat
    Omegahat <- sigres
    Omega0hat <- NULL
    loglik <- -n * r * (1 + log(2 * pi)) / 2
    for (i in 1 : p) {
      eig <- eigen(sigres[[i]])
      loglik <- loglik - ng[i] * sum(log(eig$values)) / 2
    }
    if (asy == T) {
      J <- matrix(0, p * r + p * r * (r + 1) / 2, p * r + p * r * (r + 1) / 2)
      for (i in 1 : (p - 1)) {
        for (j in 1 : (p - 1)) {
          aa <- (i - 1) * r + 1
          bb <- i * r
          cc <- (j - 1) * r + 1
          dd <- j * r
          J[aa : bb, cc : dd] <- (fracN[j] * fracN[i] 
                                 / fracN[p] * diag(r) %*% invres[[p]])
        }
        ee <- (i - 1) * r + 1
        ff <- i * r       
        J[ee : ff, ee : ff] <- (fracN[i] * diag(r) 
        %*% invres[[i]] + fracN[i] ^ 2 / fracN[p] * diag(r) %*% invres[[p]])
      }
      for (i in 1 : p) {
        aaa <- r * (p - 1) + 1 + (i - 1) * r * (r + 1) / 2
        bbb <- r * (p - 1) + i * r * (r + 1) / 2
        J[aaa : bbb, aaa : bbb] <- (0.5 * fracN[i] * 
          crossprod(expan(r), kronecker(invres[[i]], invres[[i]])) %*% expan(r))
      }
      ggg <- r + 1
      hhh <- p * r + p * r * (r + 1) / 2
      iii <- (p - 1) * r + p * r * (r + 1) / 2
      J[ggg : hhh, ggg : hhh] <- J[1 : iii, 1 : iii]
      J[1 : r, ] <- 0
      J[ggg : hhh, 1 : r] <- 0
      for (i in 1 : p) {
        J[1 : r, 1 : r] <- J[1 : r, 1 : r] + fracN[i] * diag(r) %*% invres[[i]] 
      }
      for (i in 1 : (p - 1)) {
        kkk <- i * r + 1
        lll <- (i + 1) * r
        J[1 : r, kkk : lll] <- fracN[i] * (invres[[p]] - invres[[i]])
        J[kkk : lll, 1 : r] <- J[1 : r, kkk : lll] 
      }
      temp <- ginv(J) #chol2inv(chol(J))
      temp1 <- kronecker(matrix(1, 1, (p - 1)), diag(r))
      ccc <- r + 1
      ddd <- r * p
      vargroup <- temp1 %*% tcrossprod(temp[ccc : ddd, ccc : ddd], temp1)
      covMatrix <- matrix(0, r * (p + 1), r * (p + 1))
      covMatrix[1 : ddd, 1 : ddd] <- temp[1 : ddd, 1 : ddd]
      eee <- r * p + 1
      fff <- r * (p + 1)
      covMatrix[eee : fff, eee : fff] <- vargroup
      for (i in 1 : (p - 1)) {
        mmm <- i * r + 1
        nnn <- (i + 1) * r
        covMatrix[eee : fff, mmm : nnn] <- (- temp1 %*% temp[ccc : ddd, mmm : nnn])
      }
      covMatrix[ccc : ddd, eee : fff] <- t(covMatrix[eee : fff, ccc : ddd])
      covMatrix[eee : fff, 1 : r] <- (- temp1 %*% temp[ccc : ddd, 1 : r])
      covMatrix[1 : r, eee : fff] <- t(covMatrix[eee : fff, 1 : r])
      asyFm <- matrix(sqrt(diag(covMatrix[ccc : fff, ccc : fff])), r, p)
      asySE <- asyFm
      ratio <- matrix(1, r, p) 
    }
    if (fit == T) {
      Yfit <- matrix(0, n, r)
      for (i in 1:p) {
        if (i > 1) {
          Yfit[ind[(ncum[i - 1] + 1):ncum[i]], ] <- matrix(1, ng[i], 1) %*% t(mughat[ , i])
        }
        else {
          Yfit[ind[1:ncum[1]], ] <- matrix(1, ng[1], 1) %*% t(mughat[ , 1])
        }
      }
    }
  }
  else {
    muhat <- as.matrix(mY)
    Omega0hat <- t(Gamma0hat) %*% sigY %*% Gamma0hat
    etahat = betahat = etm <- list(length = p)
    for (i in 1:p) {
      etm[[i]] <- mYg[ , i] - mY
      etahat[[i]] <- crossprod(Gammahat, etm[[i]])
      betahat[[i]] <- Gammahat %*% etahat[[i]]
    }
    Sigmahat = Omegahat = mughat = invsig <- list(length = p)
    for (i in 1:p) {
      Omegahat[[i]] <- crossprod(Gammahat, sigres[[i]]) %*% Gammahat
      Sigmahat[[i]] <- Gammahat %*% tcrossprod(Omegahat[[i]], 
                       Gammahat) + Gamma0hat %*% tcrossprod(Omega0hat, Gamma0hat)
      mughat[[i]] <- muhat + betahat[[i]]
      invsig[[i]] <- chol2inv(chol(Sigmahat[[i]]))
    }
    beta <- c()
    for (i in 1 : p) {
      beta <- cbind(beta, betahat[[i]])
    }
    betahat <- beta
    mug <- c()
    for (i in 1 : p) {
      mug <- cbind(mug, mughat[[i]])
    }
    mughat <- mug
    e11 <- crossprod(Gammahat, invsigY) %*% Gammahat
    eig2 <- eigen(e11)
    r1 <- sum(log(eig2$values))
    loglik <- - n * r * (1 + log(2 * pi)) /2 - n * logDetSigY / 2 - n * r1 / 2
    for (i in 1:p) {
      e22 <- crossprod(Gammahat, sigres[[i]]) %*% Gammahat
      eig <- eigen(e22)
      loglik <- loglik - ng[i] * sum(log(eig$values))/2
    }
    if (asy == T) {
      J <- matrix(0, p * r + p * r * (r + 1) / 2, p * r + p * r * (r + 1) / 2)
      for (i in 1 : (p - 1)) {
        for (j in 1 : (p - 1)) {
          aa <- (i - 1) * r + 1
          bb <- i * r
          cc <- (j - 1) * r + 1
          dd <- j * r
          J[aa : bb, cc : dd] <- (fracN[j] * fracN[i] 
                                           / fracN[p] * diag(r) %*% invsig[[p]])
        }
        ee <- (i - 1) * r + 1
        ff <- i * r       
        J[ee : ff, ee : ff] <- (fracN[i] * diag(r) 
        %*% invsig[[i]] + fracN[i] ^ 2 / fracN[p] * diag(r) %*% invsig[[p]])
      }
      for (i in 1 : p) {
        aaa <- r * (p - 1) + 1 + (i - 1) * r * (r + 1) / 2
        bbb <- r * (p - 1) + i * r * (r + 1) / 2
        J[aaa : bbb, aaa : bbb] <- (0.5 * fracN[i] * 
        crossprod(expan(r), kronecker(invsig[[i]], invsig[[i]])) %*% expan(r))
      }
      ggg <- r + 1
      hhh <- p * r + p * r * (r + 1) / 2
      iii <- (p - 1) * r + p * r * (r + 1) / 2
      J[ggg : hhh, ggg : hhh] <- J[1 : iii, 1 : iii]
      J[1 : r, ] <- 0
      J[ggg : hhh, 1 : r] <- 0
      for (i in 1 : p) {
        J[1 : r, 1 : r] <- J[1 : r, 1 : r] + fracN[i] * diag(r) %*% invsig[[i]] 
      }
      for (i in 1 : (p - 1)) {
        kkk <- i * r + 1
        lll <- (i + 1) * r
        J[1 : r, kkk : lll] <- fracN[i] * (invsig[[p]] - invsig[[i]])
        J[kkk : lll, 1 : r] <- J[1 : r, kkk : lll] 
      }
   
      J1 <- matrix(0, p * r + p * r * (r + 1) / 2, p * r + p * r * (r + 1) / 2)
      for (i in 1 : (p - 1)) {
        for (j in 1 : (p - 1)) {
          aa <- (i - 1) * r + 1
          bb <- i * r
          cc <- (j - 1) * r + 1
          dd <- j * r
          J1[aa : bb, cc : dd] <- (fracN[j] * fracN[i] 
                                  / fracN[p] * diag(r) %*% invres[[p]])
        }
        ee <- (i - 1) * r + 1
        ff <- i * r       
        J1[ee : ff, ee : ff] <- (fracN[i] * diag(r) 
         %*% invres[[i]] + fracN[i] ^ 2 / fracN[p] * diag(r) %*% invres[[p]])
      }
      for (i in 1 : p) {
        aaa <- r * (p - 1) + 1 + (i - 1) * r * (r + 1) / 2
        bbb <- r * (p - 1) + i * r * (r + 1) / 2
        J1[aaa : bbb, aaa : bbb] <- (0.5 * fracN[i] * 
        crossprod(expan(r), kronecker(invres[[i]], invres[[i]])) %*% expan(r))
      }
      #ggg <- r + 1
      #hhh <- p * r + p * r * (r + 1) / 2
      #iii <- (p - 1) * r + p * r * (r + 1) / 2
      J1[ggg : hhh, ggg : hhh] <- J[1 : iii, 1 : iii]
      J1[1 : r, ] <- 0
      J1[ggg : hhh, 1 : r] <- 0
      for (i in 1 : p) {
        J1[1 : r, 1 : r] <- J[1 : r, 1 : r] + fracN[i] * diag(r) %*% invres[[i]] 
      }
      for (i in 1 : (p - 1)) {
        kkk <- i * r + 1
        lll <- (i + 1) * r
        J1[1 : r, kkk : lll] <- fracN[i] * (invres[[p]] - invres[[i]])
        J1[kkk : lll, 1 : r] <- J1[1 : r, kkk : lll] 
      }
      temp1 <- ginv(J1) #chol2inv(chol(J1))
      
      
      asyFm <- matrix(sqrt(diag(temp1[1 : r * p, 1 : r * p])), r, p)
      aaaa <- p * r + p * r * (r + 1) / 2
      bbbb <- r + u * (r + p - 1 - u) + p * u * (u + 1) / 2 + (r - u) * (r - u + 1) / 2
      H <- matrix(0, aaaa, bbbb)
      for (i in 1 : (p - 1)) {
        cccc <- (i - 1) * r + 1
        dddd <- i * r
        eeee <- (i - 1) * u + 1
        ffff <- i * u
        gggg <- (p - 1) * u + 1
        hhhh <- u * (r + p - 1 - u)
        H[cccc : dddd, eeee : ffff] <- Gammahat
        H[cccc : dddd, gggg : hhhh] <- kronecker(t(etahat[[i]]), Gamma0hat)
      }
      for (i in 1 : p) {
        iiii <- r * (p - 1) + (i - 1) * r * (r + 1) / 2 + 1
        jjjj <- r * (p - 1) + i * r * (r + 1) / 2
        kkkk <- (p - 1) * u + 1
        llll <- u * (r + p - 1 - u)
        mmmm <- u * (r + p - 1 - u) + (i - 1) * u * (u + 1) / 2 + 1
        nnnn <- u * (r + p - 1 - u) + i * u * (u + 1) / 2
        oooo <- u * (r + p - 1 - u) + p * u * (u + 1) / 2 + 1 
        pppp <- u * (r + p - 1 - u) + p * u * (u + 1) / 2 + (r - u) * (r - u + 1) / 2
        H[iiii : jjjj, kkkk : llll] <- 2 * contr(r) %*% (
            kronecker(Gammahat %*% Omegahat[[i]], Gamma0hat) -
            kronecker(Gammahat, Gamma0hat %*% Omega0hat))
        H[iiii : jjjj, mmmm : nnnn] <- contr(r) %*% kronecker(Gammahat, Gammahat) %*% expan(u)
        H[iiii : jjjj, oooo : pppp] <- contr(r) %*% kronecker(Gamma0hat,
                                          Gamma0hat) %*% expan(r - u)
      }
      qqqq <- r + 1
      rrrr <- p * r + p * r * (r + 1) / 2
      ssss <- r + u * (r + p - 1 - u) + p * u * (u + 1) / 2  + (r - u) * (r - u + 1) / 2
      tttt <- (p - 1) * r + p * r * (r + 1) / 2
      uuuu <- u * (r + p - 1 - u) + p * u * (u + 1) / 2 + (r - u) * (r - u + 1) / 2
      wwww <- r * p
      xxxx <- r * p + 1
      yyyy <- r * (p + 1)
      H[qqqq : rrrr, qqqq : ssss] <- H[1 : tttt, 1 : uuuu]
      H[1 : r, qqqq : ssss] <- 0
      H[qqqq : rrrr, 1 : r] <- 0
      H[1 : r, 1 : r] <- diag(r)
      temp3 <- crossprod(H, J) %*% H
      invtemp3 <- chol2inv(chol(temp3))
      temp <- H %*% tcrossprod(invtemp3, H)
      temp2 <- kronecker(matrix(1, 1, (p - 1)), diag(r))
      vargroup <- temp2 %*% tcrossprod(temp[qqqq : wwww, qqqq : wwww], temp2)
      covMatrix <- matrix(0, yyyy, yyyy)
      covMatrix[1 : wwww, 1 : wwww] <- temp[1 : wwww, 1 : wwww]
      covMatrix[xxxx : yyyy, xxxx : yyyy] <- vargroup
      for (i in 1 : (p - 1)) {
        zzzz <- i * r + 1
        zzzzz <- (i + 1) * r
        covMatrix[xxxx : yyyy, zzzz : zzzzz] <- (- temp2 %*% temp[qqqq : wwww, zzzz : zzzzz])
      }
      covMatrix[qqqq : wwww, xxxx : yyyy] <- t(covMatrix[xxxx : yyyy, qqqq : wwww])
      covMatrix[xxxx : yyyy, 1 : r] <- (- temp2 %*% temp[qqqq : wwww, 1 : r])
      covMatrix[1 : r, xxxx : yyyy] <- t(covMatrix[xxxx : yyyy, 1 : r])
      asySE <- matrix(sqrt(diag(covMatrix[qqqq : yyyy, qqqq : yyyy])), r, p)
      ratio <- asyFm / asySE
    }
    if (fit == T) {
      Yfit <- matrix(rep(0, n * r), ncol = r)
      for (i in 1:p) {
        if (i > 1) {
          Yfit[ind[(ncum[i - 1] + 1):ncum[i]], ] <- matrix(1, ng[i], 1) %*% t(mughat[ , i])
        }
        else {
          Yfit[ind[1:ncum[1]], ] <- matrix(1, ng[1], 1) %*% t(mughat[ ,1])
        }
      }
    }
  }
  return(list(Gamma = Gammahat, Gamma0 = Gamma0hat, beta = betahat, 
              Sigma = Sigmahat, eta = etahat, Omega = Omegahat, 
              Omega0 = Omega0hat, mu = muhat, mug = mughat, loglik = loglik, 
              n = n, ng = ng, Yfit = Yfit, covMatrix = covMatrix, 
              asySE = asySE, ratio = ratio, groupInd = groupind))
}
