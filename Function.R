
library("MASS")
library("Matrix")
library("maxLik")
library("evmix")
library("np")
library("foreach")
library("doParallel")

NWEst <- function(index, y, bw, kernel_type){
  
  dist_mat <- abs(outer(index, index, "-"))
  # registerDoParallel(cores = 2) 
  # len <- length(y)
  # kernel_mat <- tryCatch({
  #   foreach(i = 1:len, .combine = rbind) %dopar% {
  #     kernel_type(dist_mat[i, ]/bw)
  #   }
  # })
  # kernel_mat <- t(kernel_mat)
  kernel_mat <- kernel_type(dist_mat/bw)
  diag(kernel_mat) <- 0
  row_sums <- rowSums(kernel_mat)
  if (any(row_sums == 0)) {
    row_sums[row_sums == 0] <- 1e-6
  }
  return(y %*% kernel_mat / row_sums)
}

ASBRM <- function(x, y, v, tau){
  
  logLik <- function(pars){
    if (ncol(x) == 2) {
      index <- x[, 1] + x[, 2] * pars[1]
    }else{
      index <- as.vector(x[, 1] + x[, -1] %*% pars[1:(ncol(x) - 1)])}
    ConExp <- NWEst(index, y, pars[ncol(x)], kernel_type = kdgaussian)
    ConExp <- pmax(pmin(ConExp, 1 - 1e-6), 1e-6)
    
    sum(c(tau * log(ConExp[y == 1]), (1 - tau) * log(1 - ConExp[y == 0])))
  }
  
  ineqA <- matrix(c(rep(0, ncol(x) - 1), 1), nrow = 1, ncol = ncol(x))
  start_value_beta <- coef(RcppArmadillo::fastLm(x, y))[-1] / coef(RcppArmadillo::fastLm(x, y))[1]
  start_value <- c( start_value_beta, sd(v)*1.06*length(x)^(-1/5) )
  
  summary(maxLik(logLik, start = start_value, method = "BFGS", constraints = list(ineqA = ineqA, ineqB = 0)))
}

ASBRM_loop <- function(x, tau, nloop) {
  
  result_est <- matrix(0, nrow = nloop, ncol = ncol(x) + 1)
  result_std <- matrix(0, nrow = nloop, ncol = ncol(x))
  
  for (i in 1:nloop) {
    
    source("DGP.r")
    res <- ASBRM(x = X, y = Y, v = V, tau)
    result_est[i,] <- c(1, res$estimate[,1])
    result_std[i,] <- res$estimate[,2]
    
    if (any(is.nan(result_std[i,]))) {next}
    # Logit <- glm(y ~ x1 + x2, family = binomial(link = "logistic"))
  }
  
  output <- list(est = result_est,
                 std = result_std)
  return(output)
}

ASBRM_contrastloop <- function(x, y, tau, nloop) {
  
  result_est_as <- matrix(0, nrow = nloop, ncol = ncol(x) + 1)
  result_std_as <- matrix(0, nrow = nloop, ncol = ncol(x))
  
  result_est_ks <- matrix(0, nrow = nloop, ncol = ncol(x) + 1)
  result_std_ks <- matrix(0, nrow = nloop, ncol = ncol(x))
  
  for (i in 1:nloop) {
    
    source("DGP.r")
    res <- ASBRM(x = X, y = Y, v = V, tau)
    result_est_as[i,] <- c(1, res$estimate[,1])
    result_std_as[i,] <- res$estimate[,2]
    if (any(is.nan(result_std[i,]))) {next}

    bw <- npindexbw(xdat = X, ydat = Y, method = "kleinspady")
    res_ks <- npindex(bws = bw, gradients = TRUE)
    result_est_ks[i,] <- c(1, res_ks$beta[2])
    result_std_ks[i,] <- sqrt(res_ks$betavcov[ncol(x), ncol(x)])
  }
  
  output <- list(est_as = result_est_as,
                 std_as = result_std_as,
                 est_ks = result_est_ks,
                 std_ks = result_std_ks)
  return(output)
}
