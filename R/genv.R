genv <- function(X, Y, Z, u, asy = TRUE, fit = TRUE, init = NULL){
  
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
    skip <-0
  }
  
  if(skip == 0){
    
    ncumx <- c()
    for (i in 1 : p) {
      ncumx[i] <- length(which(ZZ == as.numeric(levels(ZZ)[i])))
    }
    ncum <- cumsum(ncumx)
    ng <- diff(c(0, ncum))
    sortz <- sort(Z, index.return = T)
    Zs <- sortz$x
    ind <- sortz$ix
    Ys <- Y[ind, ]
    Xs <- X[ind, ]
    mY <- apply(Y, 2, mean)
    zz <- list(length = p)
    for(i in 1:p){
      if(i > 1) {
        zz[[i]] <- stats::cov(Xs[((ncum[i - 1] + 1):ncum[i]), ])
      }else{
        zz[[i]] <- stats::cov(Xs[(1:ncum[i]), ])
      }
    }
    sigres = sigYcg = sigresc <- list(length = p)
    mYg <- matrix(rep(0, r*p), ncol = p)
    mXg <- matrix(rep(0, c*p), ncol = p)
    for (i in 1 : p) {
      if(i > 1) {
        yy <- Ys[(ncum[i-1] + 1):ncum[i], ]
        xx <- Xs[(ncum[i-1] + 1):ncum[i], ]
        sigres[[i]] <- stats::cov(yy) * (ng[i] - 1)/ng[i]
        mYg[ , i] <- apply(yy, 2, mean)
        mXg[ , i] <- apply(xx, 2, mean)
        Ycg <- scale(yy, center=T, scale=F)
        Xcg <- scale(xx, center=T, scale=F)
        QZC <- diag(ng[i]) - Xcg %*% solve(t(Xcg) %*% Xcg) %*% t(Xcg)
        sigYcg[[i]] <- t(Ycg) %*% Ycg
        sigresc[[i]] <- t(Ycg) %*% QZC %*% Ycg/ng[i]
      } else {
        yy <- Ys[1:ncum[i], ]
        xx <- Xs[1:ncum[i], ]
        sigres[[i]]<- stats::cov(yy) * (ng[i] - 1)/ng[i]
        mYg[ , i] <- apply(yy, 2, mean)
        mXg[ , i] <- apply(xx, 2, mean)
        Ycg <- scale(yy, center=T, scale=F)
        Xcg <- scale(xx, center=T, scale=F)
        QZC <- diag(ng[i]) - Xcg %*% solve(t(Xcg) %*% Xcg) %*% t(Xcg)
        sigYcg[[i]] <- t(Ycg) %*% Ycg
        sigresc[[i]] <- t(Ycg) %*% QZC %*% Ycg/ng[i]
      }
    }
    sigY <- stats::cov(Y) * (n - 1)/n
    eigtemY <- eigen(sigY)$values
    logDetSigY <- log(prod(eigtemY[eigtemY > 0]))
    invsigY <- solve(sigY)
    
    sigYc <- matrix(rep(0, r*r), ncol = r)
    for (i in 1 : p){
      sigYc <- sigYc + sigYcg[[i]]/n
    }
    eigtemYc <- eigen(sigYc)$values
    logDetSigYc <- log(prod(eigtemYc[eigtemYc > 0]))
    invsigYc <- solve(sigYc)
    
    U = M <- list(length = p)
    for (i in 1 : p){
      M[[i]] <- sigresc[[i]]
      U[[i]] <-  sigYc - M[[i]]
    }
    MU <- sigYc
    tmp <- genvMU(M, U, MU, u, n, ng, p)
    
    if (!is.null(init)) {
        if (nrow(init) != r || ncol(init) != u) stop("The initial value should have r rows and u columns.")
        tmp0 <- qr.Q(qr(init), complete = TRUE)
        tmp$Gammahat <- as.matrix(tmp0[, 1:u])
        tmp$Gamma0hat <- as.matrix(tmp0[, (u+1):r])
    }
    
    Gammahat <- tmp$Gammahat
    Gamma0hat <- tmp$Gamma0hat
    covMatrix <- NULL
    asySE <- NULL
    ratio <- NULL
    Yfit <- NULL
    
    if(u == 0){
      
      Sigmahat <- sigYc
      etahat <- NULL
      Omegahat <- NULL
      Omega0hat <- sigYc
      betahat <- list(length=p)
      for (i in 1:p){
        betahat[[i]] <- matrix(rep(0, r*c), ncol=c)
      }
      muhat <- mYg
      bb <- rep(0, p)
      for (i in 1:p){
        if(i > 1){
          yy <- Ys[(ncum[i - 1]+1):ncum[i],]
          Ycg <- scale(yy, center = T, scale = F)
          bb[i] <- sum(diag(Ycg %*% solve(sigYc) %*% t(Ycg)))
        } else {
          yy <- Ys[1 : ncum[i], ]
          Ycg <- scale(yy, center = T, scale = F)
          bb[i] <- sum(diag(Ycg %*% solve(sigYc) %*% t(Ycg)))
        }
      }
      b <- sum(bb)
      loglik <- - n * r * log(2 * pi)/2 - n * logDetSigYc/2 - b/2
      if (asy == T) {
        ratio <- list(length = p)
        for(i in 1:p) {
          ratio[[i]] <- matrix(1, r, c)
        }
      }
      if (fit == T) {
        Yfit <- matrix(rep(0, n * r), ncol = r)
      }
      
    } else if (u == r) {
      
      Sigmahat <- sigresc 
      aa <- rep(0, p)
      etahat <- list(length = p)
      for (i in 1 : p){
        if(i > 1){
          yy <- Ys[(ncum[i - 1]+1):ncum[i], ]
          xx <- Xs[(ncum[i - 1]+1):ncum[i], ]
          Ycg <- scale(yy, center = T, scale = F)
          Xcg <- scale(xx, center = T, scale = F)
          etahat[[i]] <- t(Ycg) %*% Xcg %*% solve(t(Xcg) %*% Xcg)
          QZC <- diag(ng[i]) - Xcg %*% solve(t(Xcg) %*% Xcg) %*% t(Xcg)
          aa[i] <- sum(diag(QZC %*% Ycg %*% 
                              solve(sigresc[[i]]) %*% t(Ycg) %*% QZC))
        } else {
          yy <- Ys[1:ncum[i], ]
          xx <- Xs[1:ncum[i], ]
          Ycg <- scale(yy, center = T, scale = F)
          Xcg <- scale(xx, center = T, scale = F)
          etahat[[i]] <- t(Ycg) %*% Xcg %*% solve(t(Xcg) %*% Xcg)
          QZC <- diag(ng[i]) - Xcg %*% solve(t(Xcg) %*% Xcg) %*% t(Xcg)
          aa[i] <- sum(diag(QZC %*% Ycg %*% solve(sigresc[[i]])
                            %*% t(Ycg) %*% QZC))
        }
      }
      a <- sum(aa)
      Omegahat <- sigresc
      Omega0hat <- NULL
      betahat <- etahat
      muhat <- mYg
      for (i in 1:p){
        muhat[ , i] <- muhat[ , i] - betahat[[i]] %*% mXg[ , i]
      }
      loglik <- - n * r * log(2 * pi)/2 - a/2
      for (i in 1 : p){
        eig <- eigen(sigresc[[i]])
        loglik <- loglik - ng[i]*sum(log(eig$values))/2
      }
      if (asy == T) {
        asybeta = sebeta = asySEbeta = ratio <- list(length = p)
        for(i in 1:p){
          asybeta[[i]] <- n * kronecker(solve(zz[[i]]), sigresc[[i]]) / ng[i]
          sebeta[[i]] <- diag(asybeta[[i]])
          asySEbeta[[i]] <- matrix(sqrt(sebeta[[i]]), ncol = c) 
        }
        covMatrix <- asybeta
        asySE <- asySEbeta
        for(i in 1:p) {
        ratio[[i]] <- matrix(1, r, c)
        }
      }
      if (fit == T) {
        Yfit <- matrix(rep(0, n * r), ncol = r)
        for (i in 1 : p) {
          if(i > 1) {
            Yfit[ind[(ncum[i - 1] + 1):ncum[i]], ] <- rep(1, ng[i]) %*% t(muhat[ ,i]) + X[ind[(ncum[i - 1] + 1):ncum[i]],] %*% t(betahat[[i]]) 
          } else {
            Yfit[ind[1:ncum[1]], ] <- rep(1, ng[1]) %*% t(muhat[ ,1]) + X[ind[1:ncum[1]], ] %*% t(betahat[[1]]) 
          }
        }
      }
      
    } else {
      Omega0hat <- t(Gamma0hat) %*% sigYc %*% Gamma0hat
      Sigmahat = Omegahat <- list(length = p)
      for(i in 1 : p){
        Omegahat[[i]] <- t(Gammahat) %*% sigresc[[i]] %*% Gammahat
        Sigmahat[[i]] <- Gammahat %*% Omegahat[[i]] %*% t(Gammahat) + Gamma0hat %*% Omega0hat %*% t(Gamma0hat)
      }
      aa = bb <- rep(0, p)
      etahat = betahat <- list(length = p)
      for (i in 1 : p){
        if(i > 1){
          yy <- Ys[(ncum[i - 1] + 1):ncum[i], ]
          xx <- Xs[(ncum[i - 1] + 1):ncum[i], ]
          Ycg <- scale(yy, center = T, scale = F)
          Xcg <- scale(xx, center = T, scale = F)
          etahat[[i]] <- t(Gammahat) %*% t(Ycg) %*% Xcg %*% solve(t(Xcg) %*% Xcg)
          betahat[[i]] <- Gammahat %*% etahat[[i]]
          QZC <- diag(ng[i]) - Xcg %*% solve(t(Xcg) %*% Xcg) %*% t(Xcg)
          aa[i] <- sum(diag(QZC %*% Ycg %*% Gammahat %*% solve(t(Gammahat) %*% sigresc[[i]] %*% 
                                                                 Gammahat) %*% t(Gammahat) %*% t(Ycg) %*% QZC))
          bb[i] <- sum(diag(Ycg %*% Gamma0hat %*% solve(t(Gamma0hat) %*% 
                                                          sigYc %*% Gamma0hat) %*% t(Gamma0hat) %*% t(Ycg)))
        }else{
          yy <- Ys[1 : ncum[i], ]
          xx <- Xs[1 : ncum[i], ]
          Ycg <- scale(yy, center = T, scale = F)
          Xcg <- scale(xx, center = T, scale = F)
          etahat[[i]] <- t(Gammahat) %*% t(Ycg) %*% Xcg %*% solve(t(Xcg) %*% Xcg)
          betahat[[i]] <- Gammahat %*% etahat[[i]]
          QZC <- diag(ng[i]) - Xcg %*% solve(t(Xcg) %*% Xcg) %*% t(Xcg)
          aa[i] <- sum(diag(QZC %*% Ycg %*% Gammahat %*% 
                              solve(t(Gammahat) %*% sigresc[[i]] %*% Gammahat) %*% t(Gammahat) %*% t(Ycg) %*% QZC))
          bb[i] <- sum(diag(Ycg %*% Gamma0hat %*% solve(t(Gamma0hat) %*% 
                                                          sigYc %*% Gamma0hat) %*% t(Gamma0hat) %*% t(Ycg)))
        }
      }
      muhat <- mYg
      for (i in 1 : p){
        muhat[ , i]<- muhat[ , i] - betahat[[i]] %*% mXg[ , i]
      }
      a <- sum(aa)
      b <- sum(bb)
      eig2 <- eigen(t(Gamma0hat) %*% sigYc %*% Gamma0hat)
      r1 <- sum(log(eig2$values))
      loglik <- -n * r * log(2 * pi)/2 - a/2 -b/2 - n * r1/2
      for (i in 1:p){
        eig <- eigen(t(Gammahat) %*% sigresc[[i]] %*% Gammahat)
        loglik <- loglik - ng[i] * sum(log(eig$values))/2
      }
      if (asy == T) {
        asybeta = sebeta = asyFmbeta <- list(length = p)
        for(i in 1:p){
          asybeta[[i]] <- n * kronecker(solve(zz[[i]]), sigresc[[i]]) / ng[i]
          sebeta[[i]] <- diag(asybeta[[i]])
          asyFmbeta[[i]] <- matrix(sqrt(sebeta[[i]]), ncol = c) 
        }
        covMatrix <- asybeta
        asyFm <- asyFmbeta
        
        asybeta <- list(length = p)
        sebeta = asySEbeta = ratio <- list(length = p)
        bb <- matrix(rep(0, u * (r - u) * u * (r - u)), ncol = u * (r - u))
        for(i in 1 : p) {
          b <- kronecker(etahat[[i]] %*% zz[[i]] %*% t(etahat[[i]]), 
                         solve(Omega0hat)) + kronecker(Omegahat[[i]], 
                         solve(Omega0hat)) + kronecker(solve(Omegahat[[i]]), 
                         Omega0hat) - 2 * kronecker(diag(u), diag(r - u))
          bb <- bb + ng[i] * b/n
        }
        for(i in 1 : p) {
          aux <- kronecker(t(etahat[[i]]), 
                 Gamma0hat) %*% solve(bb) %*% kronecker(etahat[[i]], 
                 t(Gamma0hat))
          asybeta[[i]] <- n * kronecker(solve(zz[[i]]), Gammahat %*% 
                    Omegahat[[i]] %*% t(Gammahat)) / ng[i] + aux
          sebeta[[i]] <- diag(asybeta[[i]])
          asySEbeta[[i]] <- matrix(sqrt(sebeta[[i]]), ncol = c) 
        }
        covMatrix <- asybeta
        asySE <- asySEbeta
        for(i in 1:p){
        ratio[[i]] <- asyFm[[i]] / asySE[[i]]
        }
      }
      if (fit == T) {
        Yfit <- matrix(rep(0, n*r), ncol = r)
        for (i in 1 : p) {
          if(i > 1) {
            Yfit[ind[(ncum[i - 1] + 1):ncum[i]], ] <- rep(1, ng[i]) %*% t(muhat[ ,i]) + X[ind[(ncum[i - 1] + 1):ncum[i]], ] %*% t(betahat[[i]]) 
          } else {
            Yfit[ind[1:ncum[1]], ] <- rep(1, ng[1]) %*% t(muhat[ ,1]) + X[ind[1:ncum[1]], ] %*% t(betahat[[1]]) 
          }
        }
      }
      
    }
    
    return(list(Gamma = Gammahat, Gamma0 = Gamma0hat,
                beta = betahat, Sigma = Sigmahat,
                eta = etahat, Omega = Omegahat, Omega0 = Omega0hat, 
                mu = muhat, loglik = loglik, n = n, ng = ng, Yfit = Yfit,
                covMatrix = covMatrix, asySE = asySE, ratio = ratio, groupInd = groupind))
  }
}