reps <- 1000 # total number of simulations # 1000
upper_clamp <- 3
lower_clamp <- 0
n <- 1000 # sample size
R_synthetic <- 200 # number of synthetic samples
ep <- 1 # privacy guarantee parameter
alpha <- .05 # confidence level

R <- 3
stdmin <- 10^-8
stdmax <- 10

population_mu <- 1 # population mean
population_sigma <- 1 # population sd


# 8 * max(1/(ep/3)*log(R/stdmin/(alpha/4)), 1/(ep/3)*log(log(stdmax/stdmin)/log(2)/(alpha/4))) 
# = 381.4789
### This means that we need n to be larger than 381.4789



###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
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
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

# show progress bar
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)



###################################################
maxi <- function(v) {
  l <- 1
  for(i in 2:length(v)) {
    if(v[l] < v[i]) {
      l <- i
    }
  }
  return(l)
}

variance <- function(x, m) {
  return((1/(length(x) - 1))*sum((x - m)^2))
}


pub_interval <- function(db, a) {
  m <- mean(db)
  radius <- sqrt(variance(db, m)/length(db))*qnorm(1 - a/2)
  return(c(m - radius, m + radius))
}

pub_range <- function(db) {
  return(c(min(db), max(db)))
}

pub_histogram_learner <- function(db, bins) {
  #outputs a normalized histogram of db separated by the intervals in bins
  db <- sort(db)
  probs <- rep(0, length(bins) - 1)
  db_i <- 1
  while(db[db_i] < bins[1]) {
    # stop(paste("wrong! smaller values exist", min(db), bins[1]))
    db_i <- db_i + 1
    if(db_i > length(db)) {
      return(probs/sum(probs))
    }
  }
  for(i in 1:length(probs)) {
    while(db[db_i] < bins[i + 1]) {
      probs[i] <- probs[i] + 1
      db_i <- db_i + 1
      if(db_i > length(db)) {
        return(probs/sum(probs))
      }
    }
  }
  # stop(paste("wrong! larger values exist", max(db), bins[length(bins)]))
  return(probs/sum(probs))
}

priv_histogram_learner <- function(db, bins, e) {
  probs <- pub_histogram_learner(db, bins)
  return(probs + rlaplace(length(probs), 0, 2/e/length(db)))
}

priv_std <- function(db, a, e, stdmin, stdmax) {
  bins_base <- floor(log2(stdmin) - 2)
  bins <- 2^(bins_base:ceiling(log2(stdmax) + 2))
  y <- 1:floor(length(db)/2)
  for(i in 1:length(y)) {
    y[i] <- abs(db[2*i] - db[2*i - 1])
  }
  
  l <- maxi(priv_histogram_learner(y, bins, e))
  return(2^((l + bins_base - 1) + 2))
}

priv_mean <- function(db, e, std, r) {
  rn <- ceiling(r/std)
  bins_base <- -rn
  bins <- ((bins_base:(rn + 1)) - .5)*std
  
  l <- maxi(priv_histogram_learner(db, bins, e))
  return((l + bins_base - 1)*std)
}

priv_range <- function(db, a1, a2, e1, e2, stdmin, stdmax, r) {
  priv_std_val <- priv_std(db, a1, e1, stdmin, stdmax)
  priv_mean_val <- priv_mean(db, e2, priv_std_val, r)
  
  radius <- 4*priv_std_val*sqrt(log(length(db)/a2))
  return(c(priv_mean_val - radius, priv_mean_val + radius))
}


priv_vadhan <- function(db, a0, a1, a2, a3, e1, e2, e3, stdmin, stdmax, r) { #stdmin, stdmax, r) {
  n <- length(db)
  #rpiv_range is the most significant source of error
  xrange <- priv_range(db, a3/2, a3/2, e3/2, e3/2, stdmin, stdmax, r)
  xmin <- xrange[1]
  xmax <- xrange[2]
  xdist <- xmax - xmin
  
  #clamp
  db[db < xmin] <- xmin
  db[db > xmax] <- xmax
  
  mean_var <- xdist/(e1*n)
  priv_mean <- mean(db) + rlaplace(1, 0, mean_var)
  if(priv_mean < xmin) {
    priv_mean = xmin
  } else if(priv_mean > xmax) {
    priv_mean = xmax
  }
  
  var_var <- xdist^2/(e2*(n - 1))
  #priv_var <- public_var + extra_var + lap_noise
  priv_var <- variance(db, priv_mean) + var_var*log(1/a2) + rlaplace(1, 0, var_var)
  if(priv_var < 0 || priv_var > stdmax) {
    priv_var = stdmax
  } 
  
  #mean_var*log(1/a1) is the second most significant source of error
  priv_radius <- sqrt(priv_var/n)*qt(1 - a0/2, n - 1) + mean_var*log(1/a1)
  
  return(c(priv_mean - priv_radius, priv_mean + priv_radius))
}

priv_vadhan_ci <- function(db, a, e, stdmin, stdmax, R) {
  n <- length(db)
  if (n < 8 * max(1/(e/3)*log(R/stdmin/(a/4)), 1/(e/3)*log(log(stdmax/stdmin)/log(2)/(a/4)))){
    return(c(-R, R))
  }
  return(priv_vadhan(db, a/4, a/4, a/4, a/4, e/3, e/3, e/3, stdmin, stdmax, R))
}




###################################################
results <- foreach(
  rep = 1:reps, 
  .combine = 'cbind',
  .packages=c('rmutil'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data_randomness <- rnorm(n)
  db <- population_sigma * data_randomness + population_mu
  CI <- priv_vadhan_ci(db, alpha, ep, stdmin, stdmax, R)
  current_covered <- (CI[1] <= population_mu & CI[2] >= population_mu)
  current_width <- CI[2] - CI[1]
  rbind(current_covered, current_width, CI[1], CI[2])
}

covered <- results[1,]
width <- results[2,]
result_stat <- data.frame("coverage"=mean(covered),
                          "se_coverage"=sqrt(mean(covered) * (1 - mean(covered)) / reps),
                          "mean_width"=mean(width),
                          "se_width"=sqrt(var(width) / reps),
                          "median_width"=median(width))
result_stat
#  coverage se_coverage mean_width   se_width median_width
#     0.994  0.00244213   3.313589 0.01277874     3.412893


