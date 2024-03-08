###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "tidyverse",
  "purrr"
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

reps <- 1000
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
###################################################
source("implementation.R")


pub_test <- function(data){
  wilcox.test(data[data[, 2] == 0, 1], data[data[, 2] == 1, 1])$p.value
}

### grid search. (Note that this is only for possibly optimal result, 
###               and in real-world experiment we cannot do this.)

# here we manually tune the best parameter for each settings.
n1 <- 50
n2 <- 50
m_list <- c(3,4,5,6,7,8)
alpha0_list <- c(0.3,0.2,0.1,0.05)
m_list <- c(4)
alpha0_list <- c(0.2)
for(m in m_list){
  for(alpha0 in alpha0_list){
    p.values <- foreach(
      rep = 1:reps,
      .combine = 'rbind',
      .packages = c('purrr')
      # .options.snow=opts
    ) %dopar% {
      set.seed(rep)
      data <- data.frame(X = cbind(c(runif(n1), rbeta(n2, shape1=2, shape2=5)), c(rep(0, n1), rep(1, n2))))
      test.of.tests(data, pub_test, epsilon = 1, m = m, alpha0 = alpha0)$p.value
    }
    print(c(m, alpha0, mean(p.values < 0.05)))
  }
}
sim_TOT_power <- mean(p.values < 0.05)
# (20,80)
# [1] 3.000 0.300 0.239
# [1] 3.000 0.200 0.281
# [1] 3.000 0.100 0.332
# [1] 3.000 0.050 0.293
# [1] 4.00 0.30 0.28
# [1] 4.000 0.200 0.315
# [1] 4.000 0.100 0.328
# [1] 4.000 0.050 0.306
# [1] 5.00 0.30 0.27
# [1] 5.000 0.200 0.325
# [1] 5.000 0.100 0.351 *
# [1] 5.000 0.050 0.285
# [1] 6.000 0.300 0.265
# [1] 6.000 0.200 0.302
# [1] 6.000 0.100 0.285
# [1] 6.000 0.050 0.268
# [1] 7.000 0.300 0.278
# [1] 7.00 0.20 0.31
# [1] 7.000 0.100 0.248
# [1] 7.000 0.050 0.273
# [1] 8.000 0.300 0.229
# [1] 8.000 0.200 0.284
# [1] 8.000 0.100 0.196
# [1] 8.000 0.050 0.203

# (30, 70)
# [1] 3.000 0.300 0.278
# [1] 3.000 0.200 0.341
# [1] 3.000 0.100 0.369
# [1] 3.000 0.050 0.342
# [1] 4.000 0.300 0.344
# [1] 4.000 0.200 0.376
# [1] 4.000 0.100 0.411 *
# [1] 4.000 0.050 0.345
# [1] 5.000 0.300 0.329
# [1] 5.000 0.200 0.366
# [1] 5.000 0.100 0.399
# [1] 5.000 0.050 0.337
# [1] 6.000 0.300 0.343
# [1] 6.0 0.2 0.4
# [1] 6.000 0.100 0.366
# [1] 6.000 0.050 0.333
# [1] 7.000 0.300 0.329
# [1] 7.000 0.200 0.371
# [1] 7.000 0.100 0.347
# [1] 7.000 0.050 0.314
# [1] 8.000 0.300 0.302
# [1] 8.000 0.200 0.336
# [1] 8.000 0.100 0.318
# [1] 8.000 0.050 0.323

# (50, 50)
# [1] 3.000 0.300 0.357
# [1] 3.000 0.200 0.419
# [1] 3.000 0.100 0.457
# [1] 3.000 0.050 0.429
# [1] 4.000 0.300 0.406
# [1] 4.000 0.200 0.475 *
# [1] 4.000 0.100 0.473
# [1] 4.000 0.050 0.413
# [1] 5.000 0.300 0.427
# [1] 5.000 0.200 0.449
# [1] 5.000 0.100 0.463
# [1] 5.000 0.050 0.377
# [1] 6.000 0.300 0.391
# [1] 6.000 0.200 0.451
# [1] 6.000 0.100 0.428
# [1] 6.000 0.050 0.371
# [1] 7.000 0.300 0.394
# [1] 7.000 0.200 0.405
# [1] 7.000 0.100 0.383
# [1] 7.00 0.05 0.31
# [1] 8.000 0.300 0.356
# [1] 8.00 0.20 0.37
# [1] 8.00 0.10 0.35
# [1] 8.000 0.050 0.317

p.values <- foreach(
  rep = 1:reps,
  .combine = 'rbind',
  .packages = c('purrr'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data <- data.frame(X = cbind(c(runif(n1), rbeta(n2, shape1=2, shape2=5)), c(rep(0, n1), rep(1, n2))))
  pub_test(data)
}
sim_nonDP_power <- mean(p.values < 0.05)

p.values <- foreach(
  rep = 1:reps,
  .combine = 'rbind',
  .packages = c('purrr'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data <- data.frame(X = cbind(c(runif(n1), runif(n2)), c(rep(0, n1), rep(1, n2))))
  pub_test(data)
}
sim_nonDP_typeI <- mean(p.values < 0.05)


p.values <- foreach(
  rep = 1:reps,
  .combine = 'rbind',
  .packages = c('purrr'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data <- data.frame(X = cbind(c(runif(n1), runif(n2)), c(rep(0, n1), rep(1, n2))))
  test.of.tests(data, pub_test, epsilon = 1, m = m, alpha0 = alpha0)$p.value
}
sim_TOT_typeI <- mean(p.values < 0.05)

c(sim_nonDP_power, sim_nonDP_typeI, sim_TOT_power, sim_TOT_typeI)
# 0.821 0.060 0.351 0.042   #(20, 80), (5, 0.1)
# 0.902 0.055 0.411 0.043   #(30, 70), (4, 0.1)
# 0.958 0.058 0.475 0.035   #(50, 50), (5, 0.2)
