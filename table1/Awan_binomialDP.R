# ## The latest development version from Github
# install.packages("devtools")
# devtools::install_github("tranntran/binomialDP")
library(binomialDP)

set.seed(1001)
#set up the parameters
n = 100                      #number of trials
truth = 0.2                 #true theta
reps = 1000                #sample size
ep = 1                      #epsilon
de = 0.0                   #delta
b = exp(-ep)                #b
q = 2*de*b/(1-b+2*de*b)     #q
alpha = 0.05                #level alpha

sdp_fun = function(seed, N, ep, theta, n){
  B = qbinom(seed, size = n, prob = theta)
  s1 = (B) + (1/ep)*N[1]
  return(s1)
}

CITwoUpper = rep(NA, 1000)
CITwoLower = rep(NA, 1000)
for (rep in 1:1000){
  set.seed(rep)
  seed = runif(1, 0, 1)
  N = rtulap(n = 1, m = 0, b = b, q = q)
  Z = sdp_fun(seed, N, ep, theta_true, n_true)
  CITwoLower[rep] = CITwoSide(alpha = 0.05, Z, size = n, b = b, q = q)[1]
  CITwoUpper[rep] = CITwoSide(alpha = 0.05, Z, size = n, b = b, q = q)[2]
}
width = CITwoUpper - CITwoLower
covered = (CITwoUpper >= truth) & (CITwoLower <= truth)

print('Asymptotically Unbiased Two-sided Confidence Intervals (coverage)')
print(c(mean(covered), sqrt(mean(covered)*(1-mean(covered))/reps)))
print('Asymptotically Unbiased Two-sided Confidence Intervals (width)')
print(c(mean(width), sqrt(var(width)/reps)))

# [1] 0.947000000 0.007084561
# [1] 0.1632295427 0.0003630982
