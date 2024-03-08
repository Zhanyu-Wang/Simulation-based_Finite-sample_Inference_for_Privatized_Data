library(tidyverse)
library(purrr)


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
# Requires helper functions: dp.binom.p.val (at bottom of page)

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
  
  # # Run the Awan & Slavkovic (2018) private binomial test
  # Z <- rtulap(n = 1, m = sum(sub_tests <= alpha0), b = exp(-epsilon))
  # p_val <- 1 - dp.binom.p.val(n = m, alpha0 = alpha0, epsilon = epsilon, Z = Z)
  ### Instead, we use Gaussian DP (mu=epsilon) and the sensitivity being 1.
  Z <- rnorm(n = 1, mean = sum(sub_tests <= alpha0), sd = 1 / epsilon)
  p_val <- 1 - dp.binom.p.val(n = m, alpha0 = alpha0, sd = 1 / epsilon, Z = Z)
  return(list("Z" = Z, "p.value" = p_val))
}








######################################################
#### Helper Functions
######################################################


# Compute p-value given test statistic, Z
# Algorithm 1 in Awan & Slavkovic (2018)
dp.binom.p.val <- function(n, alpha0, sd, Z){
  F_underbar <- 0:n %>%
    map_dbl(~ pnorm(Z, mean=.x, sd=sd))
  B_underbar <- 0:n %>%
    map_dbl(~ choose(n, .x) * alpha0^.x * (1 - alpha0)^(n - .x))
  return(sum(F_underbar * B_underbar))
}

