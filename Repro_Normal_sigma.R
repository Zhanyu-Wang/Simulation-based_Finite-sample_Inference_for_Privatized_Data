# simulation settings
reps <- 1000 # total number of simulations # 1000
upper_clamp <- 3
lower_clamp <- 0
n <- 100 # sample size
R_synthetic <- 200 # number of synthetic samples
ep <- 1 # privacy guarantee parameter
alpha <- .05 # confidence level

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

scoreMCE <- function(optim_par, sa, data_randomness, privacy_noises, dp_statistic){
  mu <- optim_par
  t <- (floor(alpha * (R_synthetic + 1)) + 1)
  synth <- sdp_vec(data_randomness, privacy_noises, sa, mu)
  synth <- rbind(synth, dp_statistic)
  D_synth <- depth.Mahalanobis(synth, synth)
  D_sort <- sort(D_synth)
  # return(D_sort[t] - D_synth[R_synthetic + 1]) 
  r = sum(D_sort<=D_synth[R_synthetic+1])
  s = r + D_synth[R_synthetic+1]
  return(-s)
}

scoreMCE2 <- function(optim_par, data_randomness, privacy_noises, dp_statistic){
  mu <- optim_par[1]
  sa <- optim_par[2]
  return(scoreMCE(mu, sa, data_randomness, privacy_noises, dp_statistic))
}

accept <- function(proposed_sa, data_randomness, privacy_noises, n, dp_statistic) {
  ### optim fails with larger n...
  opt <-
    optim(
      par = 1,
      fn = scoreMCE,
      lower = 0,
      upper = 10,
      method = "L-BFGS-B",
      sa = proposed_sa,
      data_randomness = data_randomness,
      privacy_noises = privacy_noises,
      dp_statistic = dp_statistic,
    )
  # mu is proposed and we find its corresponding sa to try to 'hide' the original statistic in the synthetic statistics.
  return ((-opt$value) >= (floor((alpha)*(R_synthetic+1))+1))
}


getCI <- function(dp_statistic, data_randomness, privacy_noises, n) {
  # we try to find a theta_m in the confidence interval
  theta_m = 1
  if (!accept(theta_m, data_randomness, privacy_noises, n, dp_statistic)){
    opt = optim(
      par = c(1, 1),
      fn = scoreMCE2,
      lower = c(-5, 0),
      upper = c(5, 5),
      method = "L-BFGS-B",
      data_randomness = data_randomness,
      privacy_noises = privacy_noises,
      dp_statistic = dp_statistic,
    )
    if(opt$value <= 0)
      theta_m <- opt$par[2]
    else {
      print("failed to find a point in CI")
      break
    }
  }
  ### binary search:
  theta_0 <- 0
  theta_1 <- 10
  tol <- 10 ^ -4
  
  ##get lower bound:
  lower <- theta_0
  middle <- theta_m
  while (middle - lower > tol) {
    proposal <- (middle + lower) / 2
    if (accept(proposal, data_randomness, privacy_noises, n, dp_statistic))
      middle <- proposal
    else{
      lower <- proposal
    }
  }
  
  ##get upper bound:
  middle <- theta_m
  upper <- theta_1
  while (upper - middle > tol) {
    proposal <- (upper + middle) / 2
    if (accept(proposal, data_randomness, privacy_noises, n, dp_statistic))
      middle <- proposal
    else{
      upper <- proposal
    }
  }
  
  return(c(lower, upper))
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
  rbind(current_covered, current_width)
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
