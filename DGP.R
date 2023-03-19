
library("MASS")

DGP <- function(n, p, beta, data_type){
  
  if (data_type == 1) {
    corx <- array(0.1, dim = c(p, p))
    diag(corx) <- 1
    Xdat <- mvrnorm(n, rep(0, p), corx)
  }else if (data_type == 2 && p == 2) {
    X1 <- rnorm(n)
    X2 <- (X1 + 2*rnorm(n))/sqrt(5) + 1
    Xdat <- cbind(X1, X2)
  }else{
    X1 <- rnorm(n)
    X2 <- (X1 + 2*rnorm(n))/sqrt(5) + 1
    X3 <- rnorm(n)^2/sqrt(2)
    Xdat <- cbind(X1, X2, X3)
  }
  
  Vdat <- as.vector(Xdat %*% beta)
  Edat <- rnorm(n)
  Ydat <- ifelse( Vdat + Edat > 0, 1, 0)
  
  dat_all = list(Y = Ydat, X = Xdat, V = Vdat)
  
  return(dat_all)
}

data <- DGP(N, P, beta, data_type)
Y <- data$Y
V <- data$V
X <- data$X
