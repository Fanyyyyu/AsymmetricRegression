
source("Function.R")

N <- 300 ##The number of samples (100, 300 and 1000)
P <- 2
beta <- c(2, -1)
data_type <- 1 ## "1" represents balanced data, "2" represents unbalanced data
nloop <- 200 

set.seed(12)
tau_values <- c(0.2, 0.5, 0.8)
Result_list <- list()
ResultSE_list <- list()

for (tau in tau_values) {
  Result <- ASBRM_loop(X, tau, nloop)
  Result_list[[as.character(tau)]] <- Result$est[, 2]
  ResultSE_list[[as.character(tau)]] <- Result$std[, 1]
}

mean_lst <- list()
for (i in seq_along(Result_list)) {
  mean_lst[[i]] <- mean(Result_list[[i]])
}
mean_lst

# filename <- paste0("N=", N, ",", "data_type=", data_type, ".png")
# png(filename)
# boxplot(Result_list, names = tau_values, xlab = expression(tau), ylab = expression(beta), cex.axis = 1, cex.lab = 1.2, adj = 0.9)
# dev.off()



