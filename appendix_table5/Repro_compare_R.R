n <- 100
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
  "rmutil",
  "tictoc"
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

sizes <- c(50, 50) # n1, n2

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
  # N = rnorm(n = 2, mean = 0, sd = 1)
  ep_m <- ep * sqrt(ep_m_ratio)
  ep_u <- ep * sqrt(1 - ep_m_ratio)
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
  # N    = matrix(rnorm(n = R * 2, mean = 0, sd = 1), ncol = 2)
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

  ep_m <- ep * sqrt(ep_m_ratio)
  ep_u <- ep * sqrt(1 - ep_m_ratio)
  v <- 1 / ep_m ^ 2 # var of rnorm(1, 0, 1 / ep_m)
  D_synth = (u_tilde - m_tilde * (n - m_tilde) / 2 - v / 2) /
    sqrt(m_tilde2 * (n - m_tilde2) * (n + 1) / 12  + 2 * (n / ep_u)^2 +
           (n - 2 * m_tilde2)^2 * v / 4 + v^2 / 2) # v^2/2 = E[Z^4]/4 + v^2/4 - E[Z^2]v/2
  s = rank(D_synth, ties.method="max")[R+1]
  return(s)
}


####################################     
# search m
accept = function(seed, N, ep, ep_m_ratio, DP_n1_est, DP_test_stat){
  ep_m <- ep * sqrt(ep_m_ratio)
  DP_n1_est <- s_dp[1] 
  DP_test_stat <- s_dp[2] 
  n1_low = min(floor(n / 2), max(1, floor(DP_n1_est - qnorm(1 - alpha0, 0, 1 / ep_m))))
  n1_high = min(floor(n / 2), max(1, ceiling(DP_n1_est + qnorm(1 - alpha0, 0, 1 / ep_m))))
  n1_cand = seq(n1_low, n1_high)
  R <- length(seed[,1])

  for(n1 in n1_cand){
    val = scoreMCE(seed, N, ep, ep_m_ratio, n1, s_dp)
    if(val >= (floor((alpha)*(R+1))+1)){
      return(n1)
    }
  }
  return(NA)
}

# search m
for(R in c(100,200,500,1000,2000,5000)){
  for(small_n in c(20)){
    sizes <- c(small_n, 100 - small_n)
    print(sizes)
    for(ep_m_ratio in c(0.2)){ 
      tic()
      typeI_error <- foreach(
        i = 1:reps,
        .combine = 'rbind',
        .packages = c('rmutil', 'ddalpha'),
        .options.snow=opts
      ) %dopar% {
        set.seed(i)
        n1 <- min(sizes)
        
        seed <- runif(n, 0, 1)
        N <- rnorm(n = 2, mean = 0, sd = 1)
        s_dp <- sdp_fun(seed, N, ep, ep_m_ratio, n1)
        
        seed <- matrix(runif(R * n, 0, 1), nrow = R, ncol = n)
        N <- matrix(rnorm(R * 2, mean = 0, sd = 1), nrow = R, ncol = 2)
        accept(seed, N, ep, ep_m_ratio, s_dp[1], s_dp[2])
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
        N <- rnorm(n = 2, mean = 0, sd = 1)
        s_dp <- DP_mw_test(data, ep, ep_m_ratio, N)
        
        seed <- matrix(runif(R * n, 0, 1), nrow = R, ncol = n)
        N <- matrix(rnorm(R * 2, mean = 0, sd = 1), nrow = R, ncol = 2)
        accept(seed, N, ep, ep_m_ratio, s_dp[1], s_dp[2])
      }
      time_elapsed = toc()
      print(c(ep_m_ratio, mean(is.na(typeI_error)), mean(is.na(simu_power)), R, time_elapsed$toc - time_elapsed$tic))
      
      # write.csv(cbind(typeI_error, simu_power), file=paste("search_m_test_R=", R, ".csv", sep=''))
    }
  }
}
   
# search m
# 0.200   0.023  0.397 100.000   8.227
# 0.200   0.030  0.441 200.000  12.849
# 0.200   0.035  0.464 500.000  31.645
# 0.200   0.039  0.477 1000.000  63.392
# 0.200   0.042  0.488 2000.000  128.635
# 0.200   0.042  0.494 5000.000  321.166



####################################     
# fix m
accept = function(seed, N, ep, ep_m_ratio, DP_n1_est, DP_test_stat){
  ep_m <- ep * sqrt(ep_m_ratio)
  DP_n1_est <- s_dp[1] 
  DP_test_stat <- s_dp[2] 
  R <- length(seed[,1])
  
  n1_cand <- c(20)
  for(n1 in n1_cand){
    val = scoreMCE(seed, N, ep, ep_m_ratio, n1, s_dp)
    if(val >= (floor((alpha)*(R+1))+1)){
      return(n1)
    }
  }
  return(NA)
}

for(R in c(100,200,500,1000,2000,5000)){
  for(small_n in c(20)){
    sizes <- c(small_n, 100 - small_n)
    print(sizes)
    for(ep_m_ratio in c(0.2)){ 
      tic()
      typeI_error <- foreach(
        i = 1:reps,
        .combine = 'rbind',
        .packages = c('rmutil', 'ddalpha'),
        .options.snow=opts
      ) %dopar% {
        set.seed(i)
        n1 <- min(sizes)
        
        seed <- runif(n, 0, 1)
        N <- rnorm(n = 2, mean = 0, sd = 1)
        s_dp <- sdp_fun(seed, N, ep, ep_m_ratio, n1)
        
        seed <- matrix(runif(R * n, 0, 1), nrow = R, ncol = n)
        N <- matrix(rnorm(R * 2, mean = 0, sd = 1), nrow = R, ncol = 2)
        accept(seed, N, ep, ep_m_ratio, s_dp[1], s_dp[2])
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
        N <- rnorm(n = 2, mean = 0, sd = 1)
        s_dp <- DP_mw_test(data, ep, ep_m_ratio, N)
        
        seed <- matrix(runif(R * n, 0, 1), nrow = R, ncol = n)
        N <- matrix(rnorm(R * 2, mean = 0, sd = 1), nrow = R, ncol = 2)
        accept(seed, N, ep, ep_m_ratio, s_dp[1], s_dp[2])
      }
      time_elapsed = toc()
      print(c(ep_m_ratio, mean(is.na(typeI_error)), mean(is.na(simu_power)), R, time_elapsed$toc - time_elapsed$tic))
      
      # write.csv(cbind(typeI_error, simu_power), file=paste("fix_m_test_R=", R, ".csv", sep=''))
    }
  }
}

# fix m
# 0.200   0.052   0.519 100.000   3.195 
# 0.200   0.054   0.526 200.000   2.270 
# 0.200   0.052   0.530 500.000   3.645 
# 0.200    0.049    0.529 1000.000    4.511 
#  0.200    0.050    0.527 2000.000    8.219 
#  0.200    0.052    0.532 5000.000   21.554 
