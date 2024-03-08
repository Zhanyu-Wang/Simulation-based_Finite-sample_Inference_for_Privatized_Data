sizes <- c(50, 50) # n1, n2
n <- sum(sizes)
alpha0 <- 0.000
alpha <- 0.05 - alpha0
ep <- 1
reps = 1e3
R <- 1000

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

DP_mw_test <- function(data, ep, ep_m_ratio, N){
  # N = rlaplace(n = 2, m = 0, s = 1)
  ep_m <- ep * ep_m_ratio
  ep_u <- ep * (1 - ep_m_ratio)
  X <- data[data[, 2] == 0, 1] 
  Y <- data[data[, 2] == 1, 1]
  n1 <- length(X)
  n2 <- length(Y)
  m <- min(n1, n2)
  m_tilde <- m + N[1] * 1 / ep_m
  u_tilde <- mw(X, Y) + N[2] * n / ep_u
  return(c(m_tilde, u_tilde))
}

sdp_fun = function(seed, N, ep, ep_m_ratio, n1){
  # seed = runif(n)
  n2 <- n - n1
  data <- data.frame(X = cbind(seed, c(rep(0, n1), rep(1, n2))))
  return(DP_mw_test(data, ep, ep_m_ratio, N))
}

sdp_vec = function(seed, N, ep, ep_m_ratio, n1){
  # seed = matrix(runif(R * n), ncol = n)
  # N    = matrix(rlaplace(n = R * 2, m = 0, s = 1), ncol = 2)
  R <- length(seed[,1])
  results <- matrix(rep(0, R*2), nrow = R)
  for(i in 1:R){
    results[i,] <- sdp_fun(seed[i,], N[i,], ep, ep_m_ratio, n1)
  }
  return(results)
}

scoreMCE = function(seed, N, ep, ep_m_ratio, n1, s_dp){
  synth = sdp_vec(seed, N, ep, ep_m_ratio, n1)
  synth = rbind(synth, s_dp)
  # D_synth = depth.Mahalanobis(synth, synth)
  m_tilde <- synth[,1]
  m_tilde2 <- pmin(pmax(m_tilde, 1), as.integer(n/2))
  u_tilde <- synth[,2]

  ep_m <- ep * ep_m_ratio
  ep_u <- ep * (1 - ep_m_ratio)
  v <- 2 / ep_m ^ 2 # var of rlaplace(1, 0, 1 / ep_m)
  D_synth = (u_tilde - m_tilde * (n - m_tilde) / 2 - v / 2) /
    sqrt(m_tilde2 * (n - m_tilde2) * (n + 1) / 12 + 2 * (n / ep_u)^2 +
           (n - 2 * m_tilde2)^2 * v / 4 - v^2 / 4 + 6 / ep_m^4 / 4)
  s = rank(D_synth, ties.method="max")[R+1]
  return(s)
}


accept = function(seed, N, ep, ep_m_ratio, DP_n1_est, DP_test_stat){
  ep_m <- ep * ep_m_ratio
  DP_n1_est <- s_dp[1] 
  DP_test_stat <- s_dp[2] 
  n1_low = min(floor(n / 2), max(1, floor(DP_n1_est - qlaplace(1 - alpha0, 0, 1 / ep_m))))
  n1_high = min(floor(n / 2), max(1, ceiling(DP_n1_est + qlaplace(1 - alpha0, 0, 1 / ep_m))))
  n1_cand = seq(n1_low, n1_high)
  
  for(n1 in n1_cand){
    val = scoreMCE(seed, N, ep, ep_m_ratio, n1, s_dp)
    if(val >= (floor((alpha)*(R+1))+1)){
      return(TRUE)
    }
  }
  return(FALSE)
}

for(m in c(20, 30, 50)){
  sizes <- c(m, 100 - m)
  print(sizes)
  for(ep_m_ratio in c(0.2, 0.3, 0.4)){ 
    typeI_error <- foreach(
      i = 1:reps,
      .combine = 'rbind',
      .packages = c('rmutil', 'ddalpha'),
      .options.snow=opts
    ) %dopar% {
      set.seed(i)
      n1 <- min(sizes)
      
      seed <- runif(n, 0, 1)
      N <- rlaplace(n = 2, m = 0, s = 1)
      s_dp <- sdp_fun(seed, N, ep, ep_m_ratio, n1)
      
      seed <- matrix(runif(R * n, 0, 1), nrow = R, ncol = n)
      N <- matrix(rlaplace(R * 2, m = 0, s = 1), nrow = R, ncol = 2)
      !accept(seed, N, ep, ep_m_ratio, s_dp[1], s_dp[2])
    }
    
    simu_power <- foreach(
      i = 1:reps,
      .combine = 'rbind',
      .packages = c('rmutil'),
      .options.snow=opts
    ) %dopar% {
      set.seed(i)
      n1 <- min(sizes)
      n2 <- max(sizes)
      data <- data.frame(X = cbind(c(runif(n1), rbeta(n2, shape1=2, shape2=5)), c(rep(0, n1), rep(1, n2))))
      N <- rlaplace(n = 2, m = 0, s = 1)
      s_dp <- DP_mw_test(data, ep, ep_m_ratio, N)
      
      seed <- matrix(runif(R * n, 0, 1), nrow = R, ncol = n)
      N <- matrix(rlaplace(R * 2, m = 0, s = 1), nrow = R, ncol = 2)
      !accept(seed, N, ep, ep_m_ratio, s_dp[1], s_dp[2])
    }
    print(c(alpha0, ep_m_ratio, mean(typeI_error), mean(simu_power)))
  }
}

####################################
# alpha = 0.05 ######
# R=1000
# sizes <- c(20, 80)  

# alpha0 ep_m  typeI power
# 0.000 0.200 0.055 0.205
# 0.000 0.300 0.048 0.239 
# 0.000 0.400 0.044 0.219

# sizes <- c(30, 70)  
# 0.000 0.200 0.040 0.363
# 0.000 0.300 0.037 0.380 
# 0.000 0.400 0.038 0.320

# sizes <- c(50, 50)  
# 0.000 0.200 0.045 0.606
# 0.000 0.300 0.043 0.547
# 0.000 0.400 0.041 0.446




