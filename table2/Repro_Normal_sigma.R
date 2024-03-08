# simulation settings
reps <- 1000 # total number of simulations # 1000
upper_clamp <- 3
lower_clamp <- 0
n <- 100 # sample size
R_synthetic <- 200 # number of synthetic samples
ep <- 1 # privacy guarantee parameter
alpha <- .05 # confidence level
tol <- 10^-4

population_mu <- 1 # population mean
population_sigma <- 1 # population sd

###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha",
  "parallelly"
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
n.cores <- parallelly::availableCores(constraints = "connections")
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

# show progress bar
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


start_time <- Sys.time()
#################################################
sdp_fun <- function(data_randomness, privacy_noises, sa, mu) {
  n <- length(data_randomness)
  data <- sa * data_randomness + mu
  data_clamp <- pmax(pmin(data, upper_clamp), lower_clamp)
  s1 <- mean(data_clamp) + (upper_clamp - lower_clamp) / (n * ep) * privacy_noises[1]
  s2 <- var(data_clamp) + (upper_clamp - lower_clamp) ^ 2 / (n * ep) * privacy_noises[2]
  return(c(s1, s2))
}

sdp_vec <- function(data_randomness, privacy_noises, sa, mu) {
  # data_randomness should be (R_synthetic x n)
  # privacy_noises should be (R_synthetic x 2)
  n <- dim(data_randomness)[2]
  data <- sa * data_randomness + mu
  data_clamp <- pmax(pmin(data, upper_clamp), lower_clamp)
  s1 <- apply(data_clamp, 1, mean) + (upper_clamp - lower_clamp) / (n * ep) * privacy_noises[, 1]
  s2 <- apply(data_clamp, 1, var) + (upper_clamp - lower_clamp) ^ 2 / (n * ep) * privacy_noises[, 2]
  return(cbind(s1, s2))
}

scoreMCE2 <- function(optim_par, data_randomness, privacy_noises, dp_statistic){
  mu <- optim_par[2]
  sa <- optim_par[1]
  synth <- sdp_vec(data_randomness, privacy_noises, sa, mu)
  synth <- rbind(synth, dp_statistic)
  D_synth <- depth.Mahalanobis(synth, synth)
  r = rank(D_synth, ties.method="max")[R_synthetic+1]
  s = r + D_synth[R_synthetic+1]
  return(-s)
}

accept <- function(optim_par, data_randomness, privacy_noises, dp_statistic, search_left_bound, search_right_bound) {
  # check whether the proposal (optim_par) is acceptable
  proposed_result <- scoreMCE2(optim_par, data_randomness, privacy_noises, dp_statistic)
  if ((-proposed_result) >= (floor((alpha)*(R_synthetic+1))+1)){
    return(optim_par)
  }

  opt <-
    optim(
      par = optim_par,
      fn = scoreMCE2,
      lower = c(search_left_bound, -10),
      upper = c(search_right_bound, 10),
      method = "L-BFGS-B",
      data_randomness = data_randomness,
      privacy_noises = privacy_noises,
      dp_statistic = dp_statistic,
    )
  if ((-opt$value) >= (floor((alpha)*(R_synthetic+1))+1)){
    return(opt$par)
  } else {
    # cannot find any param that is acceptable
    return(c(NA, NA))
  }
}


getCI <- function(dp_statistic, data_randomness, privacy_noises, n) {
  theta_lower_bound <- tol
  theta_upper_bound <- 10
  # we try to find a theta_m in the confidence interval
  theta_m = 1
  nuisance_value = 1
  optim_par <- c(theta_m, nuisance_value)
  search_left_bound <- theta_lower_bound
  search_right_bound <- theta_upper_bound
  optim_par <- accept(optim_par, data_randomness, privacy_noises, dp_statistic, search_left_bound, search_right_bound)
  if (is.na(optim_par[1])){
    print("failed to find a point in CI")
    return(c(0, 0))
  }
  theta_accepted_middle_val <- optim_par[1]
  
  ## get CI left end
  left_lower <- theta_lower_bound
  left_upper <- theta_accepted_middle_val - tol
  optim_par_new <- optim_par
  while (left_upper - left_lower > tol){
    left_middle <- (left_lower + left_upper) / 2
    # check left half
    search_left_bound <- left_lower
    search_right_bound <- left_middle
    optim_par_new[1] <- search_left_bound # (search_left_bound + search_right_bound) / 2
    optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_left_bound, search_right_bound)
    if (is.na(optim_par_new[1])){
      # try it again with different init
      optim_par_new <- optim_par
      optim_par_new[1] <- (search_left_bound + search_right_bound) / 2
      optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_left_bound, search_right_bound)

      if (is.na(optim_par_new[1])){
        # exclude the left half
        left_lower <- search_right_bound
        optim_par_new <- optim_par
      }
      else {
        # find something in the left half
        left_upper <- optim_par_new[1] - tol
      }
    }
    else{
      # find something in the left half
      left_upper <- optim_par_new[1] - tol
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
    optim_par_new[1] <- search_right_bound # (search_left_bound + search_right_bound) / 2
    optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_left_bound, search_right_bound)
    if (is.na(optim_par_new[1])){
      # try it again with different init
      optim_par_new <- optim_par
      optim_par_new[1] <- (search_left_bound + search_right_bound) / 2
      optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_left_bound, search_right_bound)

      if (is.na(optim_par_new[1])){
        # exclude the left half
        right_upper <- search_left_bound
        optim_par_new <- optim_par
      }
      else {
        # find something in the right half
        right_lower <- optim_par_new[1] + tol
      }
    }
    else{
      # find something in the right half
      right_lower <- optim_par_new[1] + tol
    }
  }
  
  return(c(left_lower, right_upper))
}





results <- foreach(
  rep = 1:reps, 
  .combine = 'cbind',
  .packages=c('ddalpha'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data_randomness <- rnorm(n)
  privacy_noises <- rnorm(n = 2)
  dp_statistic <- sdp_fun(data_randomness, privacy_noises, population_sigma, population_mu)
  
  data_randomness <- matrix(rnorm(n * R_synthetic), ncol = n, nrow = R_synthetic)
  privacy_noises <- matrix(rnorm(R_synthetic * 2), ncol = 2, nrow = R_synthetic)
  CI <- getCI(dp_statistic, data_randomness, privacy_noises, n)
  current_covered <- (CI[1] <= population_sigma & CI[2] >= population_sigma)
  current_width <- CI[2] - CI[1]
  rbind(current_covered, current_width, CI[1], CI[2])
}

description_str <- paste("repro_sa-mu=", population_mu, 
                         "_sigma=", population_sigma, 
                         "_alpha=", alpha,
                         "_ep=", ep,
                         "_R=", R_synthetic,
                         "_n=", n,
                         "_clamp_l_u=", lower_clamp,
                         "_", upper_clamp,
                         "_reps=", reps,
                         sep='')
txt_name <- paste(description_str,
                  ".results.csv", sep='')
write.table(t(results), txt_name, sep=",", row.names=FALSE, col.names=TRUE)

covered <- results[1,]
width <- results[2,]
result_stat <- data.frame("coverage"=mean(covered),
                          "se_coverage"=sqrt(mean(covered) * (1 - mean(covered)) / reps),
                          "mean_width"=mean(width),
                          "se_width"=sqrt(var(width) / reps),
                          "median_width"=median(width),
                          "Total time"=Sys.time() - start_time)

txt_name <- paste(description_str,
                  ".stat.csv", sep='')
write.table(t(result_stat), txt_name, sep=",", row.names=TRUE, col.names=TRUE)
print(result_stat)

stopCluster(cl)
