sizes <- c(20, 80)
alpha <- 0.05
ep <- 1
reps = 1e3

###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "tidyverse",
  "rmutil"
)
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

# loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i,
      character.only = TRUE
    )
  )
}

# use doSNOW for parallel computing
n.cores <- min(124, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
###################################################



mw <- function(X, Y) {
  n1 <- length(X)
  n2 <- length(Y)
  X_and_Y <- c(X, Y)
  U1 <- n1*n2 + (n1*(n1 + 1)/2) - sum(rank(X_and_Y)[1:n1], na.rm = TRUE)
  U2 <- n1*n2 + (n2*(n2 + 1)/2) - sum(rank(X_and_Y)[(n1+1):(n1+n2)], na.rm = TRUE)
  U <- min(U1, U2)
  return(U)
}

DP_mw_test <- function(data, ep, ep_m_ratio){
  ep_m <- ep * ep_m_ratio
  ep_u <- ep * (1 - ep_m_ratio)
  X <- data[data[, 2] == 0, 1] 
  Y <- data[data[, 2] == 1, 1]
  n1 <- length(X)
  n2 <- length(Y)
  m <- min(n1, n2)
  m_tilde <- m + rlaplace(n = 1, m = 0, s = 1 / ep_m)
  u_tilde <- mw(X, Y) + rlaplace(n = 1, m = 0, s = (n1 + n2) / ep_u)
  return(c(m_tilde, u_tilde))
}

calc_null_val <- function(n1, n2, ep, ep_m_ratio){
  null_stat_values <- rep(0, 1000) 
  for(i in 1:1000) {
    set.seed(i+10000000)
    data <- data.frame(X = cbind(c(runif(n1), runif(n2)), c(rep(0, n1), rep(1, n2))))
    null_stat_values[i] <- DP_mw_test(data, ep, ep_m_ratio)[2]
  }
  return(null_stat_values)
}

for(m in c(20, 30, 50)){
  sizes <- c(m, 100 - m)
  print(sizes)
  for(ep_m_ratio in c(0.2, 0.3, 0.4)){
    typeI_error <- foreach(
      i = 1:reps,
      .combine = 'rbind',
      .packages = c('rmutil'),
      .options.snow=opts
    ) %dopar% {
      set.seed(i)
      n1 <- sizes[1]
      n2 <- sizes[2]
      n <- n1 + n2
      data <- data.frame(X = cbind(c(runif(n)), c(rep(0, n1), rep(1, n2))))
      m_and_u_tilde <- DP_mw_test(data, ep, ep_m_ratio)
      m_est <- min(n, max(0, m_and_u_tilde[1]))
      m_est <- min(n - m_est, m_est)
      m_est <- ceiling(m_est)
      if(m_est == 0){
        rej <- FALSE
      } else {
        null_values <- calc_null_val(m_est, n - m_est, ep, ep_m_ratio)
        rej_threhold <- quantile(null_values, alpha)
        rej <- m_and_u_tilde[2] < rej_threhold
      }
      rej
    }
    
    simu_power <- foreach(
      i = 1:reps,
      .combine = 'rbind',
      .packages = c('rmutil'),
      .options.snow=opts
    ) %dopar% {
      set.seed(i)
      n1 <- sizes[1]
      n2 <- sizes[2]
      n <- n1 + n2
      data <- data.frame(X = cbind(c(runif(n1), rbeta(n2, shape1=2, shape2=5)), c(rep(0, n1), rep(1, n2))))
      m_and_u_tilde <- DP_mw_test(data, ep, ep_m_ratio)
      m_est <- min(n, max(0, m_and_u_tilde[1]))
      m_est <- min(n - m_est, m_est)
      m_est <- ceiling(m_est)
      if(m_est == 0){
        rej <- FALSE
      } else {
        null_values <- calc_null_val(m_est, n - m_est, ep, ep_m_ratio)
        rej_threhold <- quantile(null_values, alpha)
        rej <- m_and_u_tilde[2] < rej_threhold
      }
      rej
    }
    print(c(ep_m_ratio, mean(typeI_error), mean(simu_power)))
  }
}


####################################
# alpha = 0.05 ######
# sizes <- c(20, 80) 

# ep_m    typeI   power
# 0.200 0.112 0.437
# 0.300 0.083 0.381
# 0.400 0.069 0.325

# sizes <- c(30, 70)
# 0.200 0.065 0.522
# 0.300 0.056 0.474
# 0.400 0.053 0.388

# sizes <- c(50, 50)
# 0.200 0.051 0.641
# 0.300 0.052 0.587
# 0.400 0.055 0.510
##########################################















