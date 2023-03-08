reps <- 1000

R <- 200
n = 500
alpha=.05
ep=1
shape1=1/2
shape2=1/2
beta0=1/2
beta1=beta1_true=2

set.seed(42)
dir.create(file.path('.', 'results/logistic'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path('.', 'results/summary'), showWarnings = FALSE, recursive = TRUE)




###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "tictoc",
  "MASS"
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
n.cores <- min(63, parallel::detectCores() - 1)
cl <- makeSOCKcluster(n.cores)
registerDoSNOW(cl)

# # show progress bar
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
# opts <- list()

###################################################


expit = function(x){
  return(exp(x)/(1+exp(x)))
}


ObjPertLog = function(ep,c,X,Y){
  n = length(Y)
  t = 1/4
  c = max(c, t/(2*n)/(exp(ep)-1))
  x_dim = dim(X)[2]-1
  lambda = (1/4)*(x_dim+1)
  ep_new = ep - log(1+t/(2*n*c))
  
  sim_norm = rgamma(n=1, shape=x_dim+1, rate=ep_new/2)
  sim_vector = rnorm(n=x_dim+1, m=0, s=1)
  b = sim_vector*sim_norm/sqrt(sum(sim_vector^2))
    
  obj = function(beta){
    return((1/n)*sum(log(1+exp(-Y * X%*%beta)))   + (c/n)*sum(beta^2) + sum(b*beta)/n)
  }
  grad = function(beta){
    return( -(1/n) * t(Y * expit(-Y * X%*%beta)) %*% X    + (2*c/n)*beta + b/n)
  }
  min = optim(par = rep(0, x_dim+1), fn = obj, gr = grad, method = "BFGS")
  
  if(min$convergence!=0){
    print("Obj-Pert did not converge") 
  }
  return(c(ep_new, min$par))
}###   END OBJ  PERT


PrivSPDMat = function(original_matrix, sensitivity, ep, c){
  d = dim(original_matrix)[1]
  norm = rgamma(n=1, shape=d^2, rate=ep/sensitivity)
  vector = rnorm(n=d^2, m=0, s=1)
  eta = matrix(vector*norm/sqrt(sum(vector^2)), nrow=d)
  new_matrix = original_matrix + eta
  new_matrix = (new_matrix + t(new_matrix)) / 2
  ev <- eigen(new_matrix)
  new_ev <- pmax(2*c, ev$values)
  new_matrix2 <- ev$vectors %*% diag(new_ev) %*% t(ev$vectors)
  return(new_matrix2)
}


n_list = c(100, 200, 500, 1000, 2000)
ep_list = c(0.1, 0.3, 1, 3, 10)

for(ep in ep_list){
  for(n in n_list){
    filename2 = paste("./results/summary/ERM_Logistic", 
                      "-shape1_", shape1, "-shape2_", shape2, 
                      "-beta0_", beta0, "-beta_1", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps, ".csv", sep="")
    print(filename2)
    if (file.exists(filename2)) {
      next
    }
    CIs <- foreach(
      rep = 1:1000,
      .combine = 'rbind',
      .packages=c('ddalpha', 'tictoc'),
      .options.snow=opts
    ) %dopar% {
      tic()
      set.seed(rep)
      
      ep1 = ep/2
      ep2 = ep1/2
      ep3 = ep1/2
      
      x = rbeta(n,shape1,shape2)
      x2 = 2*x-1
      U = runif(n,min=0,max=1)
      Y = (U<=expit(beta0+beta1*x2)) * 2 - 1
      X = matrix(c(rep(1,n),x2),nrow=n,ncol=2)
      
      q = 0.85
      t = 1/4
      c = t/(2*n)/(exp(ep1 * (1-q))-1)
      
      x_dim <- dim(X)[2]-1
      private_estimate = ObjPertLog(ep=ep1,c=c,X=X,Y=Y)
      if(sum(private_estimate^2) > 200) private_estimate[2:(x_dim+2)] = rep(0, x_dim+1)
      ep1_for_simulation = ep1 * q
      private_beta = private_estimate[2:3]
      
      new_matrix = X
      multiplier = expit(Y * X%*%private_beta) * expit(-Y * X%*%private_beta)
      new_matrix[,1] = multiplier * new_matrix[,1]
      new_matrix[,2] = multiplier * new_matrix[,2]
      hessian_estimate = t(new_matrix) %*% X / n + 2 * c * diag(2)
      hessian_sensitivity = 1/(2*n)
      print(hessian_estimate)
      private_hessian = PrivSPDMat(hessian_estimate, hessian_sensitivity, ep2, c)
      
      new_matrix = X
      multiplier = expit(-Y * X%*%private_beta)^2
      new_matrix[,1] = multiplier * new_matrix[,1]
      new_matrix[,2] = multiplier * new_matrix[,2]
      cov_estimate = t(new_matrix) %*% X / n - 4 * c^2 * private_beta %*% t(private_beta)
      cov_sensitivity = 2/n * expit(sqrt(sum(private_beta^2)))^2
      private_cov = PrivSPDMat(cov_estimate, cov_sensitivity, ep3, c)
      
      m <- 10000
      sim_G <- mvrnorm(n = m, mu = c(0,0), Sigma = private_cov)
      
      sim_norm = rgamma(n=m, shape=x_dim+1, rate=ep1_for_simulation/2)
      sim_b = matrix(rnorm(n=(x_dim+1)*m, m=0, s=1), nrow=m)
      sim_norm2 =  (sim_norm / sqrt(rowSums(sim_b^2)))
      sim_b[,1] = sim_b[,1] * sim_norm2
      sim_b[,2] = sim_b[,2] * sim_norm2
      
      sim_noise = (sim_G + sim_b / sqrt(n)) %*% solve(private_hessian) / sqrt(n)
      beta0_CI = private_beta[1] + quantile(sim_noise[,1], c(alpha/2, 1-alpha/2))
      beta1_CI = private_beta[2] + quantile(sim_noise[,2], c(alpha/2, 1-alpha/2))
      beta1_CI = pmax(-10, pmin(10, beta1_CI))
      CI = as.numeric(beta1_CI)
      time_elapsed = toc()
      result = data.frame(c(CI, time_elapsed$toc - time_elapsed$tic))
      colnames(result) = "1"
      rownames(result) = c("l", "h", "t")
      t(result)
    }
    
    write.csv(CIs, filename2)
  }
}
stopCluster(cl)
