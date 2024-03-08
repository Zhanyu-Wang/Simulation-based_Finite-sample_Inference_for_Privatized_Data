library(tidyverse)
library(purrr)
library(poibin)


######################################################
#### Run The Test of Tests for Given Data Set
######################################################
#
# `data` is a dataframe with n rows (i.e., each observation is a row)
# `pub_test` is a function that takes a dataframe as input and returns the 
#     public p-value
# `sub_samp_sizes` is an optional length-m vector containing the desired sizes of
#     each of the subsamples (if NA, points are distributed as evenly as possible)
#
# Requires helper functions: rtulap, ptulap, and dp.binom.p.val (at bottom of page)

test.of.tests <- function(data, pub_test, epsilon, m, alpha0, sub_samp_sizes = NA){
  # Partition the data into random subsamples
  if(!is.na(sub_samp_sizes[1])){
    sub_samples <- sample(rep(1:m, sub_samp_sizes))
  }
  # If no subsample sizes provided, set sizes to be floor(n/m) and ceiling(n/m)
  else{
    n <- nrow(data)
    sub_samples <- sample(rep(1:m, ceiling(n/m))[1:n])
  }
  
  sub_tests <- rep(NA, m)
  for(i in 1:m){
    # Compute p-value in each subsample
    sub_tests[i] <- try(pub_test(data[sub_samples == i,]), silent = T)
    
    # If test results in an error, sample from Unif(0,1)
    if(class(sub_tests[i]) == "character"){
      sub_tests[i] <- runif(1)
    }
  }
  
  # Run the Awan & Slavkovic (2018) private binomial test
  Z <- rtulap(n = 1, m = sum(sub_tests <= alpha0), b = exp(-epsilon))
  p_val <- 1 - dp.binom.p.val(n = m, alpha0 = alpha0, epsilon = epsilon, Z = Z)
  return(list("Z" = Z, "p.value" = p_val))
}







######################################################
#### Compute the Theoretical Power of the Test of Tests
######################################################
#
# `effect` is a scalar representing the effect size of the test
# `pub_power` is a function that takes an as input effect, n, and alpha -- which
#     represent the effect size, sample size, and significance level, respectively
#     -- and returns the public test's power
# `sub_samp_sizes` is an optional length-m vector containing the desired sizes of
#     each of the subsamples (if NA, points are distributed as evenly as possible)
#
# Requires helper functions: ptulap, dp.binom.p.val, dp.binom.alt.prob (at bottom of page)

power.test.of.tests <- function(effect, pub_power, epsilon, n, m, alpha0, alpha = 0.05, 
                                sub_samp_sizes = NA, ...){
  # If no subsample sizes provided, set sizes to be floor(n/m) and ceiling(n/m)
  if(is.na(sub_samp_sizes[1])){
    sub_samp_sizes <- c(rep(ceiling(n/m), floor(n - floor(n/m)*m)),
                        rep(floor(n/m), ceiling(m - n + floor(n/m)*m)))
  }
  
  # Create vector of the power of the public test in each subsample
  thetas <- rep(NA, m)
  for(i in 1:m){
    thetas[i] <- try(pub_power(n = sub_samp_sizes[i], effect = effect, 
                               alpha = alpha0, ...), silent = T)
    # If pub_pow encounters an error or returns 0, set power to alpha0
    if(class(thetas[i]) == "character" | thetas[i] == 0){
      thetas[i] <- alpha0
    }
  }
  
  # Find critical value (1-alpha quantile of null distribution)
  tol <- 2*log(100)/epsilon
  critical_value <- optimize(function(x){abs(dp.binom.p.val(m, alpha0, epsilon, x) - (1 - alpha))},
                             interval = c(-tol, m+tol))$minimum
  
  # Return the power (upper tail probability of distribution under alt. hypothesis)
  return(1 - dp.binom.alt.prob(m, thetas, epsilon, critical_value))
}







######################################################
#### Find the m,alpha0 that Maximize Power with Known Effect Size
######################################################
#
# `effect` is a scalar representing the effect size of the test
# `pub_power` is a function that takes an as input effect, n, and alpha -- which
#     represent the effect size, sample size, and significance level, respectively
#     -- and returns the public test's power
# `m_grid` is an optional vector of m's to search over (if NA, defaults to a 
#     reasonable set; the default set may yield slow performance for large n)
#
# Requires helper functions: ptulap, dp.binom.p.val, dp.binom.alt.prob (at bottom of page)
#
# Currently only sub-sample sizes of (approximately) n/m supported


optimal.m.alpha0 <- function(effect, pub_power, epsilon, n, m_grid = NA, 
                             alpha = 0.05, ...){
  # If no grid for m provided, assign the default grid discussed in the paper
  if(is.na(m_grid[1])){
    m_grid <- c(1:sqrt(n), floor(n/rev(1:(sqrt(n)+1))))
    m_grid <- m_grid[!duplicated(m_grid)]
  }
  
  # An function to efficiently compute the power of ToT with balanced sub-samples
  efficient_power <- function(alpha0, m, effect, pub_power, epsilon, n,  
                              alpha = 0.05, sub_samp_sizes = NA, ...){
    pub_pow1 <- try(pub_power(effect = effect, n = ceiling(n/m), alpha = alpha0, ...))
    if(class(pub_pow1) == "character" | pub_pow1 == 0){ pub_pow1 <- alpha0}
    pub_pow2 <- try(pub_power(effect = effect, n = floor(n/m), alpha = alpha0, ...))
    if(class(pub_pow2) == "character" | pub_pow2 == 0){ pub_pow2 <- alpha0}
    thetas <- c(rep(pub_pow1, floor(n - floor(n/m)*m)),
                rep(pub_pow2, ceiling(m - n + floor(n/m)*m)))
    tol <- 2*log(100)/epsilon
    critical_value <- optimize(function(x){abs(dp.binom.p.val(m, alpha0, epsilon, x) - (1 - alpha))},
                               interval = c(-tol, m+tol), tol = 5e-4)$minimum
    return(1 - dp.binom.alt.prob(m, thetas, epsilon, critical_value))
  }
  
  # Function to find the alpha0 that maximizes power for a given m
  optimal_alpha0 <- function(m, n, effect, epsilon, alpha, pub_power, ...){
    opt <- optimize(f = efficient_power, interval = c(0,1), maximum = T, tol = 5e-3,
                    m = m, effect = effect, pub_power = pub_power, epsilon = epsilon,
                    n = n, alpha = alpha, ...)
    return(c(opt$objective, opt$maximum, m))
  }
  
  # For each m in grid, find the alpha0 that maximizes power
  x <- lapply(X = m_grid, FUN = optimal_alpha0, n = n, effect = effect, 
              epsilon = epsilon, alpha = alpha, pub_power = pub_power, ...)
  
  # Find the m,alpha0 combination that yields the maximum overall power and return summary
  max_x <- x[[which.max(as.numeric(sapply(x,"[[",1)))]]
  return(list("power" = max_x[1], "alpha0" = max_x[2], "m" = as.integer(max_x[3])))
}







######################################################
#### Find m,alpha0 with Unknown Effect Size
######################################################
#
# `rho` is the target power of the test
# `pub_power` is a function that takes an as input effect, n, and alpha -- which
#     represent the effect size, sample size, and significance level, respectively
#     -- and returns the public test's power
# `m_grid` is an optional vector of m's to search over (if NA, defaults to a 
#     reasonable set; the default set may yield slow performance for large n)
# `effect_grid` is an optional vector of effect sizes to search over (if NA,  
#     defaults to a reasonable set)
#
# Requires helper functions: ptulap, dp.binom.p.val, dp.binom.alt.prob (at bottom of page)
#
# Currently only sub-sample sizes of (approximately) n/m supported

practical.m.alpha0 <- function(rho, pub_power, epsilon, n, m_grid = NA, 
                               effect_grid = NA, alpha = 0.05, ...){
  # If no grid for m provided, assign the default grid discussed in the paper
  if(is.na(m_grid[1])){
    m_grid <- c(1:sqrt(n), floor(n/rev(1:(sqrt(n)+1))))
    m_grid <- m_grid[!duplicated(m_grid)]
  }
  # If no grid for effect provided, assign a default grid
  if(is.na(effect_grid[1])){
    effect_grid <- 2^c(-7:-1, seq(0,4,0.5))
  }
  
  # Perform a binary search for the smallest effect that gives power greater than rho
  i <- (length(effect_grid)+1) %/% 2; power <- 0; effect <- Inf
  while(length(effect_grid) > 1){
    opt <- optimal.m.alpha0(effect = effect_grid[i], pub_power = pub_power, 
                            epsilon = epsilon, n = n, m_grid = m_grid, 
                            alpha = alpha, ...)
    if(opt$power >= rho){
      power <- opt$power; alpha0 <- opt$alpha0; m <- opt$m; effect <- effect_grid[i]
      effect_grid <- effect_grid[1:i]
      i <- (i+1) %/% 2
      
    }
    else{
      effect_grid <- effect_grid[(i+1):length(effect_grid)]
      i <- (length(effect_grid) + 1) %/% 2
    }
  }
  if(effect != effect_grid){
    opt <- optimal.m.alpha0(effect = effect_grid, pub_power = pub_power, 
                            epsilon = epsilon, n = n, m_grid = m_grid, 
                            alpha = alpha, ...)
    power <- opt$power; alpha0 <- opt$alpha0; m <- opt$m
  }
  return(list("power" = power, "alpha0" = alpha0, "m" = m, "effect" = effect_grid))
}







######################################################
#### Helper Functions
######################################################

# Sample from Tulap(m,b)
# Algorithm 2 in Awan & Slavkovic (2018)
rtulap <- function(n, m, b) {
  vect <- rep(NA, n)
  G1 <- rgeom(n, 1 - b)
  G2 <- rgeom(n, 1 - b)
  U <- runif(n, -.5, .5)
  vect <- G1 - G2 + U + m
  return(vect)
}

# CDF of Tulap(m,b)
# Definition 4.1 in Awan & Slavkovic (2018)
ptulap <- function(m, b, x){
  if (x <= round(m)) {
    (b ^ (- round(x-m))/ (1 + b)) * (b + (x - m - round(x - m) + 1/2) * (1 - b))
  } else {
    1 - ((b ^ round(x - m)) / (1 + b)) * (b + (round(x - m) - (x - m) + 1/2) * (1 - b))  }
}

# Compute p-value given test statistic, Z
# Algorithm 1 in Awan & Slavkovic (2018)
dp.binom.p.val <- function(n, alpha0, epsilon, Z){
  F_underbar <- 0:n %>%
    map_dbl(~ ptulap(.x, exp(-epsilon), Z))
  B_underbar <- 0:n %>%
    map_dbl(~ choose(n, .x) * alpha0^.x * (1 - alpha0)^(n - .x))
  return(sum(F_underbar * B_underbar))
}

# Compute probability of statistic as or more extreme than Z under H_A
# where H_A is given by a vector of thetas
dp.binom.alt.prob <- function(n, thetas, epsilon, Z){
  F_underbar <- 0:n %>%
    map_dbl(~ ptulap(.x, exp(-epsilon), Z))
  B_underbar <- 0:n %>%
    map_dbl(.f = dpoibin, pp = thetas)
  return(sum(F_underbar * B_underbar))
}
