reps = 1000 ### for coverage simulation
n_true = 100
theta_true = 0.2
alpha = 0.05
R = 200 

ep = 1
de = 0.0                   #delta
b = exp(-ep)                #b
q = 2*de*b/(1-b+2*de*b)     #q


###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach", 
  "doSNOW", 
  "ddalpha", 
  "parallelly", 
  "geometry",
  "rmutil",
  "binomialDP"
)
new.packages <- list.of.packages[!(list.of.packages %in% 
                                     installed.packages()[, "Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep = TRUE)
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
n.cores <- parallelly::availableCores(constraints = "connections")
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

###################################################

sdp_fun = function(seed, N, ep, theta, n){
  B = qbinom(seed, size = n, prob = theta)
  s1 = (B) + (1/ep)*N[1]
  return(s1)
}

sdp_vec = function(seed, N, ep, theta, n){
  #seed should be Rxn
  #N should be Rx1
  B = qbinom(seed, size = n, prob = theta)
  s1 = B + (1/ep)*N[, 1]
  return(s1)
}

scoreMCE = function(theta, n, seed, N, s_dp, ep, pivot = FALSE){
  n = floor(n)

  synth = sdp_vec(seed, N, ep = ep, theta = theta, n = n)  
  synth = c(synth, s_dp)
  # D_synth = depth.Mahalanobis(synth, synth)
  D_synth = 1 / (1 + (synth - mean(synth))^2 / var(synth))
  r = rank(D_synth, ties.method="max")[R+1]
  s = r + D_synth[R+1]
  return(-s)
}



accept = function(theta, seed, N, s_dp, ep, pivot, search_left_bound, search_right_bound){
  n = n_true
  val = scoreMCE(theta, n, seed, N, s_dp, ep, pivot)
  if(-val >= (floor((alpha)*(R+1))+1)){
    return(theta)
  }
  
  opt <-
    optim(
      par = theta,
      fn = scoreMCE,
      lower = c(search_left_bound),
      upper = c(search_right_bound),
      method = "L-BFGS-B",
      n = n, 
      seed = seed, 
      N = N, 
      s_dp = s_dp, 
      ep = ep, 
      pivot = pivot
    )
  if(-opt$val >= (floor((alpha)*(R+1))+1)){
    return(opt$par)
  }
  return(NA)
}


getCI = function(s_dp, seed, N, ep, pivot = FALSE){
  tol <- 1e-4
  theta_lower_bound <- 0
  theta_upper_bound <- 1
  # we try to find a theta_m in the confidence interval
  theta_m <- min(max(s_dp[1], 0), n_true) / n_true

  optim_par <- theta_m
  search_left_bound <- theta_lower_bound
  search_right_bound <- theta_upper_bound
  optim_par <- accept(optim_par, seed, N, s_dp, ep, pivot, search_left_bound, search_right_bound)
  if (is.na(optim_par)){
    print("failed to find a point in CI")
    return(0)
  }
  theta_accepted_middle_val <- optim_par
  
  ## get CI left end
  left_lower <- theta_lower_bound
  left_upper <- theta_accepted_middle_val - tol
  optim_par_new <- optim_par
  while (left_upper - left_lower > tol){
    left_middle <- (left_lower + left_upper) / 2
    # check left half
    search_left_bound <- left_lower
    search_right_bound <- left_middle
    optim_par_new <- search_left_bound # (search_left_bound + search_right_bound) / 2
    optim_par_new <- accept(optim_par_new, seed, N, s_dp, ep, pivot, search_left_bound, search_right_bound)
    if (is.na(optim_par_new)){
      # try it again with different init
      optim_par_new <- (search_left_bound + search_right_bound) / 2
      optim_par_new <- accept(optim_par_new, seed, N, s_dp, ep, pivot, search_left_bound, search_right_bound)

      if (is.na(optim_par_new)){
        # exclude the left half
        left_lower <- search_right_bound
        optim_par_new <- optim_par
      }
      else {
        # find something in the left half
        left_upper <- optim_par_new - tol
      }
    }
    else{
      # find something in the left half
      left_upper <- optim_par_new - tol
    }
  }
  
  ## get CI right end
  right_lower <- theta_accepted_middle_val + tol
  right_upper <- theta_upper_bound
  optim_par_new <- optim_par
  while (right_upper - right_lower > tol) {
    right_middle <- (right_lower + right_upper) / 2
    # check right half
    search_left_bound <- right_middle
    search_right_bound <- right_upper
    optim_par_new <- search_right_bound # (search_left_bound + search_right_bound) / 2
    optim_par_new <- accept(optim_par_new, seed, N, s_dp, ep, pivot, search_left_bound, search_right_bound)
    if (is.na(optim_par_new)){
      # try it again with different init
      optim_par_new <- optim_par
      optim_par_new <- (search_left_bound + search_right_bound) / 2
      optim_par_new <- accept(optim_par_new, seed, N, s_dp, ep, pivot, search_left_bound, search_right_bound)

      if (is.na(optim_par_new)){
        # exclude the left half
        right_upper <- search_left_bound
        optim_par_new <- optim_par
      }
      else {
        # find something in the right half
        right_lower <- optim_par_new + tol
      }
    }
    else{
      # find something in the right half
      right_lower <- optim_par_new + tol
    }
  }
  
  return(c(left_lower, right_upper))
}


results <- foreach(
  rep = 1:reps,
  .combine = 'cbind',
  .packages = c('ddalpha', 'binomialDP'),
  .options.snow = opts
) %dopar% {
  # {rep=1
  set.seed(rep)
  
  seed = runif(1, 0, 1)
  N = rtulap(n = 1, m = 0, b = b, q = q)
  s_dp = sdp_fun(seed, N, ep, theta_true, n_true)
  
  seed = runif(R, 0, 1)
  N = matrix(rtulap(n = R*1, m = 0, b = b, q = q), ncol = 1, nrow = R)

  CI <- getCI(s_dp, seed, N, ep)  

  covered <- (CI[1] <= theta_true & CI[2] >= theta_true)
  c(covered, CI[2]-CI[1])
}
stopCluster(cl)
covered <- results[1, ]
width <- results[2, ]

print('Mahalanobis Depth coverage')
print(c(mean(covered), sqrt(mean(covered)*(1-mean(covered))/reps)))
print('Mahalanobis Depth width')
print(c(mean(width), sqrt(var(width)/reps)))


# [1] "Mahalanobis Depth coverage and width"
# [1] 0.949000000 0.006956939
# [1] 0.1657240819 0.0005099669
