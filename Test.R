
set.seed(245)

data_type <- 1
source("DGP.r")
len1 <- length(which(Y == 0))
len2 <- length(which(Y == 1))
c(len1, len2)

Vdat <- as.vector(X %*% beta)
Edat <- rnorm(N)
Ydat <- ifelse( Vdat + Edat > 0, 1, 0)

dist_mat <- abs(outer(V, V, "-"))
bandwidth <- 0.02
kernel_content <- dist_mat / bandwidth
kernel_mat <- kdgaussian(dist_mat / bandwidth)
diag(kernel_mat) <- 0
row_sums <- rowSums(kernel_mat)
if (any(row_sums == 0)) {
  row_sums[row_sums == 0] <- 1e-5
}
Nominator <- Y %*% kernel_mat
NWEst <- Nominator / row_sums
as.vector(NWEst)
NWEst[Y == 0 & NWEst == 1] <- 1 - 1e-5
NWEst[Y == 1 & NWEst == 0] <- 1e-5
# if (any(NWEst > 1) || any(NWEst < 0) || any(is.na(NWEst)) || any(is.infinite(NWEst))) {
#   NWEst <- replace(NWEst, NWEst > 1, 1 - 1e-5)
#   NWEst <- replace(NWEst, NWEst < 0, 1e-5)
# }
put <- c(tau * log(NWEst[Y == 1]), (1 - tau) * log(1 - NWEst[Y == 0]))
put

