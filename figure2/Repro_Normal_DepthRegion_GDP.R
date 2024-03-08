# simulation settings
value_r = 100
upper_clamp <- 3
lower_clamp <- 0
n <- 100 # sample size
R_synthetic <- 200 # number of synthetic samples
ep <- 0.5 # privacy guarantee parameter
alpha <- .05 # confidence level
tol <- 10^-8

population_mu <- 1 # population mean
population_sigma <- 1 # population sd

###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha",
  "parallelly",
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
n.cores <- parallelly::availableCores(constraints = "connections")
print(n.cores)
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

# show progress bar
pb <- txtProgressBar(max = value_r, style = 3)
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

score_sa_mu <- function(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type){
  mu <- optim_par[2]
  sa <- optim_par[1]
  synth <- sdp_vec(data_randomness, privacy_noises, sa, mu)
  synth <- rbind(synth, dp_statistic)
  if(depth_type == 'mahalanobis'){
    D_synth = depth.Mahalanobis(synth, synth)
  }
  else if(depth_type == 'halfspace'){
    D_synth = depth.halfspace(synth, synth)
  }
  else if(depth_type == 'simplicial'){
    D_synth = depth.simplicial(synth, synth)
  }
  else if(depth_type == 'spatial'){
    D_synth = depth.spatial(synth, synth)
  }
  r = rank(D_synth, ties.method="max")[R_synthetic+1]
  s = r + D_synth[R_synthetic+1]
  return(-s)
}

score_mu_sa <- function(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type){
  mu <- optim_par[1]
  sa <- optim_par[2]
  synth <- sdp_vec(data_randomness, privacy_noises, sa, mu)
  synth <- rbind(synth, dp_statistic)
  if(depth_type == 'mahalanobis'){
    D_synth = depth.Mahalanobis(synth, synth)
  }
  else if(depth_type == 'halfspace'){
    D_synth = depth.halfspace(synth, synth)
  }
  else if(depth_type == 'simplicial'){
    D_synth = depth.simplicial(synth, synth)
  }
  else if(depth_type == 'spatial'){
    D_synth = depth.spatial(synth, synth)
  }
  r = rank(D_synth, ties.method="max")[R_synthetic+1]
  s = r + D_synth[R_synthetic+1]
  return(-s)
}

accept <- function(optim_par, data_randomness, privacy_noises, dp_statistic, search_lower, search_upper, 
                   nuisance_lower, nuisance_upper, depth_type, score_func) {
  # check whether the proposal (optim_par) is acceptable
  proposed_result <- score_func(optim_par, data_randomness, privacy_noises, dp_statistic, depth_type)
  if ((-proposed_result) >= (floor((alpha)*(R_synthetic+1))+1)){
    # print(c("accept! ", optim_par, unname(proposed_result)))
    return(optim_par)
  }

  opt <-
    optim(
      par = optim_par,
      fn = score_func,
      lower = c(search_lower, nuisance_lower),
      upper = c(search_upper, nuisance_upper),
      method = "L-BFGS-B",
      data_randomness = data_randomness,
      privacy_noises = privacy_noises,
      dp_statistic = dp_statistic,
      depth_type = depth_type
    )
  # print(c("accept? ", opt$par, opt$value))
  if ((-opt$value) >= (floor((alpha)*(R_synthetic+1))+1)){
    return(opt$par)
  } else {
    # cannot find any param that is acceptable
    return(c(NA, NA))
  }
}


getConfidenceInterval <- function(optim_par, dp_statistic, data_randomness, privacy_noises, search_lower, search_upper, nuisance_lower, nuisance_upper, depth_type, score_func) {
  # we try to find a theta_m in the confidence interval
  theta_lower <- search_lower
  theta_upper <- search_upper
  optim_par <- accept(optim_par, data_randomness, privacy_noises, dp_statistic, search_lower, search_upper,
                      nuisance_lower, nuisance_upper, depth_type, score_func)
  if (is.na(optim_par[1])){
    print("failed to find a point in CI")
    return(c(NA, NA))
  }
  theta_accepted_middle_val <- optim_par[1]
  # print(c(theta_accepted_middle_val, search_lower, search_upper, nuisance_lower, nuisance_upper))
  
  ## get CI left end
  left_lower <- theta_lower
  left_upper <- theta_accepted_middle_val - tol
  optim_par_new <- optim_par
  while (left_upper - left_lower > tol){
    left_middle <- (left_lower + left_upper) / 2
    # check left half
    search_lower <- left_lower
    search_upper <- left_middle
    optim_par_new[1] <- search_lower # (search_lower + search_upper) / 2
    optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_lower, search_upper,
                            nuisance_lower, nuisance_upper, depth_type, score_func)
    # print(c("getConfidenceInterval l", optim_par_new[1], search_lower, search_upper, left_upper))
    if (is.na(optim_par_new[1])){
      # try it again with different init
      optim_par_new <- optim_par
      optim_par_new[1] <- (search_lower + search_upper) / 2
      optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_lower, search_upper,
                              nuisance_lower, nuisance_upper, depth_type, score_func)
      # print(c("getConfidenceInterval l2", optim_par_new[1], search_lower, search_upper, left_upper))

      if (is.na(optim_par_new[1])){
        # exclude the left half
        left_lower <- search_upper
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
  right_upper <- theta_upper
  optim_par_new <- optim_par
  while (right_upper - right_lower > tol) {
    right_middle <- (right_lower + right_upper) / 2
    # check right half
    search_lower <- right_middle
    search_upper <- right_upper
    optim_par_new[1] <- search_upper # (search_lower + search_upper) / 2
    optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_lower, search_upper,
                            nuisance_lower, nuisance_upper, depth_type, score_func)
    # print(c("getConfidenceInterval r", optim_par_new[1], search_lower, search_upper, right_lower))
    if (is.na(optim_par_new[1])){
      # try it again with different init
      optim_par_new <- optim_par
      optim_par_new[1] <- (search_lower + search_upper) / 2
      optim_par_new <- accept(optim_par_new, data_randomness, privacy_noises, dp_statistic, search_lower, search_upper,
                              nuisance_lower, nuisance_upper, depth_type, score_func)
      # print(c("getConfidenceInterval r2", optim_par_new[1], search_lower, search_upper, right_lower))

      if (is.na(optim_par_new[1])){
        # exclude the left half
        right_upper <- search_lower
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

################################################################
### Plotting the confidence set

dp_statistic = c(1, 0.75)  
d1_lower = -10
d1_upper = 10
d2_lower = tol
d2_upper = 10

set.seed(1)
data_randomness <- matrix(rnorm(n * R_synthetic), ncol = n, nrow = R_synthetic)
privacy_noises <- matrix(rnorm(R_synthetic * 2), ncol = 2, nrow = R_synthetic)

################################################################
for(depth_type in c('mahalanobis', 'halfspace', 'simplicial', 'spatial')){
  print(depth_type)
  ##### getConfidenceGrid
  print(c(d1_lower, d1_upper, d2_lower, d2_upper))
  optim_par <- c(1, 1)
  d1_new_range <- getConfidenceInterval(optim_par, dp_statistic, data_randomness, privacy_noises, d1_lower, d1_upper, d2_lower, d2_upper, depth_type, score_mu_sa) # CI for mu
  d2_new_range <- getConfidenceInterval(optim_par, dp_statistic, data_randomness, privacy_noises, d2_lower, d2_upper, d1_lower, d1_upper, depth_type, score_sa_mu) # CI for sa
  print(c(d1_new_range, d2_new_range))
  boundary <- foreach(
    i = 1:(value_r^2), 
    .combine = 'rbind',
    .packages=c('ddalpha'),
    .options.snow=opts
  ) %dopar% {  
    # {i = 1
    d1_idx <- (i-1) %% value_r # start from 0
    d2_idx <- as.integer((i-1) / value_r)  # start from 0
    search_lower <- d1_new_range[1] + d1_idx / value_r * (d1_new_range[2] - d1_new_range[1])
    search_upper <- d1_new_range[1] + (d1_idx + 1) / value_r * (d1_new_range[2] - d1_new_range[1])
    nuisance_lower <- d2_new_range[1] + d2_idx / value_r * (d2_new_range[2] - d2_new_range[1])
    nuisance_upper <- d2_new_range[1] + (d2_idx + 1) / value_r * (d2_new_range[2] - d2_new_range[1])
    confidence_region_part <- c(NA, NA, NA, NA)
    if(!is.na(accept(c(search_lower + search_upper, nuisance_lower + nuisance_upper) / 2, data_randomness, privacy_noises, dp_statistic, search_lower, search_upper, 
                   nuisance_lower, nuisance_upper, depth_type, score_mu_sa)[1])){
      confidence_region_part <- c(search_lower, nuisance_lower, search_upper, nuisance_upper)
    }
    confidence_region_part
  }

  boundary_valid = boundary[!is.na(boundary[,1]),]
  write.csv(boundary_valid, file=paste(depth_type, value_r, ep, "boundary.csv", sep='_'))

  area = sum((boundary_valid[,3] - boundary_valid[,1]) * (boundary_valid[,4] - boundary_valid[,2]))
  print(c(depth_type, area))
  # ep=1, 100
  # 0.354115559299602 mahalanobis 0.5514926 1.1880036 0.7841644 1.5345536
  # 0.605506006334619 halfspace   0.312500  2.125000  0.750000  1.457031
  # 0.630322806823634 simplicial  0.2265625 2.1250000 0.7500000 1.4570313
  # 0.362314512180075 spatial     0.5568847 1.1849995 0.7767796 1.5360642

  # ep=0.5, 100
  # 0.766327478059884 mahalanobis 0.3612175  1.2197266  0.5502777  1.9074708
  # 0.953061589678131 halfspace   0.48370359 1.14062501 0.00000001 3.25000001
  # 0.964212564381321 simplicial  0.47957156 1.14062501 0.00000001 3.25000001
  # 0.776741354683335 spatial     0.3676376  1.2290956  0.5339965  1.9078484

  pdf(paste(depth_type, value_r, ep, "region.pdf", sep="_"), width = 5, height = 5)  
  plot(c(0., 4), c(0, 4), type = "n", xlab = "", ylab = "", main = area)
  rect(boundary_valid[,1], 
       boundary_valid[,2], 
       boundary_valid[,3], 
       boundary_valid[,4], 
       border = NA,
       col = "#1b98e0")

  invisible(dev.off())
}
# stopCluster(cl)

################################################################



