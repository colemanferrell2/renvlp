boot.eppls <- function(X1, X2, Y, u, B) {
    
    X1 <- as.matrix(X1)
    X2 <- as.matrix(X2)
    a <- dim(Y)
    n <- a[1]
    r <- a[2]
    p1 <- ncol(X1)
    p2 <- ncol(X2)
    
    fit <- eppls(X1, X2, Y, u, asy = F)
    Yfit <- matrix(1, n, 1) %*% t(fit$muY) + X1 %*% fit$beta1 + X2 %*% fit$beta2
    res <- Y - Yfit
    
    bootenv <- function(i) {
        res.boot <- res[sample(1:n, n, replace = T), ]
        Y.boot <- Yfit + res.boot
        newfit <- eppls(X1, X2, Y.boot, u, asy = F)
        return(c(newfit$beta1, newfit$beta2))
    }
    
    bootbeta <- lapply(1:B, function(i) bootenv(i))
    bootbeta <- matrix(unlist(bootbeta), nrow = B, byrow = TRUE)
    
    bootse <- matrix(apply(bootbeta, 2, stats::sd), nrow = r)
    bootse1 <- bootse[ , 1:p1]
    bootse2 <- bootse[ , (p1+1):(p1+p2)]
    
    return(list(bootse1=t(bootse1), bootse2=t(bootse2)))
    
}
