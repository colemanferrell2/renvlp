genvMU <- function(M, U, MU, u, n, ng, L, initial = NULL){
  
  p <- L 
  dimM <- dim(M[[1]])
  dimU <- dim(U[[1]])
  r <- dimM[1]
  
  if (!is.null(initial)) {
      if (nrow(initial) != r || ncol(initial) != u) stop("The initial value should have r rows and u columns.")
  }
  
  auxinit <- function(M, U, MU, u, ng, n, p, initial = NULL){
    
    tmp.MU <- eigen(MU)
    invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
    invMU2 <- sweep(tmp.MU$vectors, MARGIN = 2, 1 / sqrt(tmp.MU$values), '*') %*% t(tmp.MU$vectors)
    
    if (!is.null(initial)) {
        init <- initial
        obj1a <- c()
        for(i in 1:p){
            eig1 <- eigen(t(init) %*% M[[i]] %*% init)
            eig2 <- eigen(t(init) %*% invMU %*% init)
            obj1a[i] <- sum(log(eig1$values))*ng[i]/n + sum(log(eig2$values))/p
        }
        obj1 <- sum(obj1a)
    } else {    
    startv <- function(a){
      out <- list(length=p)
      for (i in 1:p){
        out[[i]] <- t(a)%*% U[[i]] %*% a
      }
      sout <- Reduce("+", out)
      return(sout)
    }
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]]) 
    
    obj1a <- c()
    for(i in 1:p){
      eig1 <- eigen(t(init) %*% M[[i]] %*% init)
      eig2 <- eigen(t(init) %*% invMU %*% init)
      obj1a[i] <- sum(log(eig1$values))*ng[i]/n + sum(log(eig2$values))/p
    }
    obj1 <- sum(obj1a)
    
    startv2 <- function(a){
      out <- list(length=p)
      for (i in 1:p){
        out[[i]] <- t(a)%*% invMU2 %*% tcrossprod(U[[i]], invMU2)  %*% a
      }
      sout <- Reduce("+", out)
      return(sout)
    }
    tmp2.MU <- apply(tmp.MU$vectors, 2, startv2)
    tmp3.MU <- sort(tmp2.MU, decreasing = TRUE, index.return = TRUE)
    init.MU <- as.matrix(tmp.MU$vectors[, tmp3.MU$ix[1:u]])
    obj1a <- c()
    for(i in 1:p){
      e1 <- eigen(t(init.MU) %*% M[[i]] %*% init.MU)
      e2 <- eigen(t(init.MU) %*% invMU %*% init.MU)
      obj1a[i] <- sum(log(e1$values))*ng[i]/n + sum(log(e2$values))/p	
    }
    obj2 <- sum(obj1a)
    
    if (obj2 < obj1) {
      init <- init.MU
      obj1 <- obj2
    }
    
    Ma = Ua <- list(length=p)
    for (i in 1:p){
      Ma[[i]] <- M[[i]]*ng[i]/n
      Ua[[i]] <- U[[i]]/p
    }
    M1 <- Reduce("+", Ma)
    U1 <- Reduce("+", Ua)
    tmp.M <- eigen(M1)
    startv3 <- function(a) t(a) %*% U1 %*% a
    tmp2.M <- apply(tmp.M$vectors, 2, startv3)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)	
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])
    e1 <- eigen(t(init.M) %*% M1 %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)
    obj3 <- sum(log(e1$values)) + sum(log(e2$values))	
    if (obj3 < obj1) {
      init <- init.M
      obj1 <- obj3
    }
    
    invM2 <- sweep(tmp.M$vectors, MARGIN = 2, 1 / sqrt(tmp.M$values), '*') %*% t(tmp.M$vectors)
    midmatrix <- invM2 %*% tcrossprod(U1, invM2) 
    startv4 <- function(a) t(a) %*% midmatrix %*% a
    tmp2.M <- apply(tmp.M$vectors, 2, startv4)
    tmp3.M <- sort(tmp2.M, decreasing = TRUE, index.return = TRUE)
    init.M <- as.matrix(tmp.M$vectors[, tmp3.M$ix[1:u]])				
    e1 <- eigen(t(init.M) %*% M1 %*% init.M)
    e2 <- eigen(t(init.M) %*% invMU %*% init.M)			
    obj4 <- sum(log(e1$values)) + sum(log(e2$values))	
    if (obj4 < obj1) {
      init <- init.M
      obj1 <- obj4
    }
  }
    
    
    return(list(init = init, obj1 = obj1, invMU = invMU))
    
  }
  
  auxf1 <- function(M1, U1, u, n, ng, p, x, r){
    M <- M1
    
    t2 <- crossprod(G1init[-j, ], as.matrix(M[-j, j])) / M[j, j]

    
    GUGt2 <- g + t2
    GUG <- crossprod(G1init, (M %*% G1init)) - tcrossprod(GUGt2, GUGt2) * M[j, j]
    
    invC1 <- chol2inv(chol(GUG))
    
    tmp2 <- x + t2
    tmp3 <- x + t3
    T2 <- invC1 %*% tmp2
    T3 <- invC2 %*% tmp3
    out <- ng * log(1 + M[j, j] * crossprod(tmp2, T2))/n + log(1 + invMU[j, j]
                  * crossprod(tmp3, T3))/p
    return(out)
  }
  
  auxf2 <- function(M1, U1, t2, t3, invc1, invc2, ng, n, p, x, j){
    M <- M1
    
    tmp2 <- x + t2
    tmp3 <- x + t3
    
    T2 <- invc1 %*% tmp2	
    T3 <- invc2 %*% tmp3
    out <- ng * log(1 + M[j, j] * crossprod(tmp2, T2))/n + log(1 + invMU[j, j] 
              * crossprod(tmp3, T3))/p
    return(out)
    
  }
  
  auxg1 <- function(M1, U1, u, n, ng, p, x, r){
    M <- M1
    
    t2 <- crossprod(G1init[-j, ], as.matrix(M[-j, j])) / M[j, j]
    
    GUGt2 <- g + t2
    GUG <- crossprod(G1init, (M %*% G1init)) - tcrossprod(GUGt2, GUGt2) * M[j, j]
    
    invC1 <- chol2inv(chol(GUG))
    
    tmp2 <- x + t2
    tmp3 <- x + t3
    T2 <- invC1 %*% tmp2	
    T3 <- invC2 %*% tmp3
    out <-  2 * ng * T2 / (n *
            as.numeric(1 / M[j, j] + crossprod(tmp2, T2))) + 2 * T3 /(p * 
            as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3)))
    return(out)
  }
  
  auxg2 <- function(M1, U1, t2, t3, invc1, invc2, ng, n, p, x, j){
    M <- M1
  
    
    tmp2 <- x + t2
    tmp3 <- x + t3
    
    T2 <- invc1 %*% tmp2	
    T3 <- invc2 %*% tmp3
    out <-  2 * T2 * ng / (n * 
            as.numeric(1 / M[j, j] + crossprod(tmp2, T2))) + 2 * T3 / (p * 
            as.numeric(1 / invMU[j, j] + crossprod(tmp3, T3)))	
    return(out)
    
  }
  
  
  if(u == 0){
    Gammahat <- NULL
    Gamma0hat <- diag(r)
  }else if (u == r){
    Gammahat <- diag(r)
    Gamma0hat <- NULL
  }else if (u == r - 1){
    
    maxiter = 100
    ftol = 1e-3
    initout <- auxinit(M, U, MU, u, ng, n, p, initial)
    init <- initout$init
    obj1 <- initout$obj1
    invMU <- initout$invMU 
    
    GEidx <- GE(init)
    Ginit = G1init <- init %*% solve(init[GEidx[1 : u], ])
    j <- GEidx[r]
    
    g <- as.matrix(Ginit[j, ])
    t3 <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
    
    GVGt2 <- g + t3
    GVG <- crossprod(Ginit, (invMU %*% Ginit)) - tcrossprod(GVGt2, GVGt2) * invMU[j, j]
    
    invC2 <- chol2inv(chol(GVG))
    
    fobj <- function(x) {
      res <- -2 * log(1 + sum(x^2)) 
      for(i in 1:p){
        res <- res + auxf1(M[[i]], U[[i]], u, n, ng[i], p, x, r)
      }
      return(res)
    }
    
    gobj <- function(x) {
      res <- -4 * x %*% solve(1 + sum(x^2))
      for(i in 1:p){
        res <- res + auxg1(M[[i]], U[[i]], u, n, ng[i], p, x, r)
      }
      return(res)
    }
    
    i <- 1
    while (i < maxiter) {
      
      res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
      Ginit[j, ] <- res$par
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      obj5a <- c()
      for(i in 1:p){
        e1 <- eigen(t(Gammahat) %*% M[[i]] %*% Gammahat)
        e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)		
        obj5a[i] <- sum(log(e1$values))*ng[i]/n + sum(log(e2$values))/p
      }
      obj5 <- sum(obj5a)
      
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u + 1) : r, drop = FALSE]
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
    
  }else{
    
    maxiter <- 100
    ftol <- 1e-3
    
    initout <- auxinit(M, U, MU, u, ng, n, p, initial)
    init <- initout$init
    obj1 <- initout$obj1
    
    GEidx <- GE(init)
    Ginit <- init %*% solve(init[GEidx[1 : u], ])
    
    GUG = GVG <- list(length = p)
    for (k in 1 : p){
      MU <- M[[k]] + U[[k]]
      tmp.MU <- eigen(MU)
      invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 
                     1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
      GUG[[k]] <- crossprod(Ginit, (M[[k]] %*% Ginit))	
      GVG[[k]] <- crossprod(Ginit, (invMU %*% Ginit))	
    }
    t4 <- crossprod(Ginit[GEidx[(u + 1):r],], Ginit[GEidx[(u + 1):r], ]) + diag(u)
    i <- 1
    while (i < maxiter) {
      
      for (j in GEidx[(u+1):r]) {
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 - tcrossprod(g, g)
        invt4 <- chol2inv(chol(t4))	
        
        t2 = t3 = GUGt2 = GVGt2 = invc1 = invc2 <- list(length=p)
        for(k in 1:p){
          Maux <- M[[k]]
          MU <- M[[k]] + U[[k]]
          tmp.MU <- eigen(MU)
          invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 
                         1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
          t2[[k]] <- crossprod(Ginit[-j, ], as.matrix(Maux[-j, j])) / Maux[j, j]
          t3[[k]] <- crossprod(Ginit[-j, ], as.matrix(invMU[-j, j])) / invMU[j, j]
          
          GUGt2[[k]] <- g + t2[[k]]
          GUG[[k]] <- GUG[[k]] - tcrossprod(GUGt2[[k]], GUGt2[[k]]) * Maux[j, j]
          
          GVGt2[[k]] <- g + t3[[k]]
          GVG[[k]] <- GVG[[k]] - tcrossprod(GVGt2[[k]], GVGt2[[k]]) * invMU[j, j] 
          
          invc1[[k]] <- ginv(GUG[[k]]) #chol2inv(chol(GUG[[k]]))
          invc2[[k]] <- ginv(GVG[[k]]) #chol2inv(chol(GVG[[k]]))
        }
        fobj <- function(x) {
          res <- -2 * log(1 + x %*% invt4 %*% x)
          for(k in 1:p){
            res <- res + auxf2(M[[k]], U[[k]], t2[[k]], t3[[k]], 
                               invc1[[k]], invc2[[k]], ng[k], n, p, x, j)
          }
          return(res)
        }
        
        gobj <- function(x) {
          res <- -4	* invt4 %*% x / as.numeric(1 + x %*% invt4 %*% x)
          for(k in 1:p){
            res <- res + auxg2(M[[k]], U[[k]], t2[[k]], t3[[k]], 
                               invc1[[k]], invc2[[k]], ng[k], n, p, x, j)
          }
          return(res)
        }
        
        res <- stats::optim(Ginit[j,], fobj, gobj, method = "BFGS")
        Ginit[j, ] <- res$par
        g <- as.matrix(Ginit[j, ])
        t4 <- t4 + tcrossprod(g, g)
        
        for(k in 1:p){
          Maux <- M[[k]]
          MU <- M[[k]] + U[[k]]
          tmp.MU <- eigen(MU)
          invMU <- sweep(tmp.MU$vectors, MARGIN = 2, 
                         1 / tmp.MU$values, '*') %*% t(tmp.MU$vectors)
          GUGt2[[k]] <- g + t2[[k]]
          GUG[[k]] <- GUG[[k]] + tcrossprod(GUGt2[[k]], GUGt2[[k]]) * Maux[j, j]
          
          GVGt2[[k]] <- g + t3[[k]]
          GVG[[k]] <- GVG[[k]] + tcrossprod(GVGt2[[k]], GVGt2[[k]]) * invMU[j, j] 
        }
        
      }
      a <- qr(Ginit)
      Gammahat <- qr.Q(a)
      obj5a <- c()
      for(i in 1 : p){
        e1 <- eigen(t(Gammahat) %*% M[[i]] %*% Gammahat)
        e2 <- eigen(t(Gammahat) %*% invMU %*% Gammahat)		
        obj5a[i] <- sum(log(e1$values))*ng[i]/n + sum(log(e2$values))/p
      }
      obj5 <- sum(obj5a)
      if (abs(obj1 - obj5) < ftol * abs(obj1)) {
        break
      } else {
        obj1 <- obj5
        i <- i + 1
      }
    }
    Gamma0hat <- qr.Q(a, complete = TRUE)[, (u + 1) : r]
    Gammahat <- as.matrix(Gammahat)
    Gamma0hat <- as.matrix(Gamma0hat)
  }
  return(list(Gammahat = Gammahat, Gamma0hat = Gamma0hat))
  
}

