###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "tidyverse",
  "purrr",
  "poibin"
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
source("implementation.R")

pub_test <- function(data){
  wilcox.test(data[data[, 2] == 0, 1], data[data[, 2] == 1, 1])$p.value
}

### grid search. (Note that this is only for possibly optimal result, 
###               and in real-world experiment we cannot do this.)
m_list <- c(3,4,5,6,7,8)
alpha0_list <- c(0.3,0.2,0.1,0.05)
m_list <- c(5)
alpha0_list <- c(0.1)
n1 <- 10
n2 <- 90

# here we manually tune the best parameter for each settings.
for(m in m_list){
  for(alpha0 in alpha0_list){
    p.values <- foreach(
      rep = 1:1000,
      .combine = 'rbind',
      .packages = c('purrr', 'poibin'),
      .export = c('ptulap'),
      .options.snow=opts
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
# m.    alpha0 power
# 3.000 0.300 0.158
# 3.00  0.20  0.18
# 3.000 0.100 0.214
# 3.000 0.050 0.182

# 4.000 0.300 0.185
# 4.000 0.200 0.231
# 4.000 0.100 0.239
# 4.000 0.050 0.217

# 5.000 0.300 0.212
# 5.000 0.200 0.259   # the best
# 5.000 0.100 0.252
# 5.000 0.050 0.204

# 6.000 0.300 0.224
# 6.000 0.200 0.238
# 6.000 0.100 0.205
# 6.000 0.050 0.189

# 7.000 0.300 0.218
# 7.000 0.200 0.258
# 7.000 0.100 0.198
# 7.000 0.050 0.191

# 8.000 0.300 0.192
# 8.000 0.200 0.238
# 8.000 0.100 0.136
# 8.000 0.050 0.129

# (30,70)
# 3.000 0.300 0.213
# 3.000 0.200 0.251
# 3.000 0.100 0.268
# 3.000 0.050 0.256

# 4.00  0.30  0.26
# 4.000 0.200 0.308
# 4.000 0.100 0.299
# 4.00  0.05  0.26

# 5.000 0.300 0.263
# 5.000 0.200 0.302
# 5.000 0.100 0.305
# 5.000 0.050 0.262

# 6.000 0.300 0.295
# 6.000 0.200 0.323
# 6.000 0.100 0.294
# 6.00  0.05  0.25

# 7.000 0.300 0.277
# 7.000 0.200 0.324  # the best 
# 7.000 0.100 0.304
# 7.000 0.050 0.248

# 8.000 0.300 0.282
# 8.000 0.200 0.306
# 8.00  0.10  0.26
# 8.000 0.050 0.244

# (50, 50)
# 3.000 0.300 0.188
# 3.000 0.200 0.259
# 3.000 0.100 0.294
# 3.000 0.050 0.261

# 4.000 0.300 0.273
# 4.000 0.200 0.307
# 4.000 0.100 0.328
# 4.00 0.05 0.28

# 5.000 0.300 0.293
# 5.000 0.200 0.342  # the best 
# 5.000 0.100 0.307
# 5.000 0.050 0.256

# 6.000 0.300 0.298
# 6.000 0.200 0.331
# 6.000 0.100 0.315
# 6.000 0.050 0.263

# 7.000 0.300 0.293
# 7.000 0.200 0.305
# 7.000 0.100 0.261
# 7.00 0.05 0.21

# 8.000 0.300 0.303
# 8.000 0.200 0.294
# 8.000 0.100 0.267
# 8.000 0.050 0.218

p.values <- foreach(
  rep = 1:1000,
  .combine = 'rbind',
  .packages = c('purrr', 'poibin'),
  .export = c('ptulap'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data <- data.frame(X = cbind(c(runif(n1), rbeta(n2, shape1=2, shape2=5)), c(rep(0, n1), rep(1, n2))))
  pub_test(data)
}
sim_nonDP_power <- mean(p.values < 0.05)

p.values <- foreach(
  rep = 1:1000,
  .combine = 'rbind',
  .packages = c('purrr', 'poibin'),
  .export = c('ptulap'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data <- data.frame(X = cbind(c(runif(n1), runif(n2)), c(rep(0, n1), rep(1, n2))))
  pub_test(data)
}
sim_nonDP_typeI <- mean(p.values < 0.05)


p.values <- foreach(
  rep = 1:1000,
  .combine = 'rbind',
  .packages = c('purrr', 'poibin'),
  .export = c('ptulap'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  data <- data.frame(X = cbind(c(runif(n1), runif(n2)), c(rep(0, n1), rep(1, n2))))
  test.of.tests(data, pub_test, epsilon = 1, m = m, alpha0 = alpha0)$p.value
}
sim_TOT_typeI <- mean(p.values < 0.05)

c(sim_nonDP_power, sim_nonDP_typeI, sim_TOT_power, sim_TOT_typeI)
                                                            # the best power setting
# 0.821            0.060            0.259          0.046.   # (20, 80), (5, 0.2)
# 0.902            0.055            0.324          0.051.   # (30, 70), (7, 0.2)
# 0.958            0.058            0.342          0.047    # (50, 50), (5, 0.2)
