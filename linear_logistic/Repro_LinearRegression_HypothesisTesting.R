# n=200
Delta = 2
mu=.5
tau=1
# X ~ N(mu, tau^2)
beta = c(-.5,1)### first entry is intercept, second is slope.
sa = .5
# Y ~ N(X * beta[2] + beta[1], sa^2)
alpha=.05
ep=1
R=200
reps=1000
dir.create(file.path('.', 'results/Repro_LR_HT'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path('.', 'results/compare_clamp'), showWarnings = FALSE, recursive = TRUE)

###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha"
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

pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


start_time <- Sys.time()

clamp = function(x,a,b){
  return(pmin(pmax(x,a),b))
}

sdp_vec = function(ux,uy,N,Delta,theta,ep,n){
  ### ux should be r*n dimensional N(0,1).
  ### uy should be r*n dimensional N(0,1)
  ### N should be r*5 dimensional N(0,1)
  ### theta should be dataframe with mu, tau, beta0, beta1, sa
  beta1=theta[1]
  beta0=theta[2]
  mu=theta[3]
  tau=theta[4]
  sa=theta[5]
  
  x = ux*tau+mu
  y = beta0+beta1*x+uy*sa
  xc = clamp(x,-Delta,Delta)
  yc = clamp(y,-Delta,Delta)
  xyc = clamp(x*y,-Delta^2,Delta^2)
  
  ep = ep/sqrt(5)### mu-GDP 
  xbar = apply((xc),1,mean)+2*Delta/(n*ep)*N[,1]
  ybar = apply((yc),1,mean)+2*Delta/(n*ep)*N[,2]
  x2bar = apply((xc^2),1,mean)+Delta^2/(n*ep)*N[,3]
  y2bar = apply((yc^2),1,mean)+Delta^2/(n*ep)*N[,4]
  xybar = apply((xyc),1,mean)+2*Delta^2/(n*ep)*N[,5]
  
  return(data.frame(xbar,ybar,x2bar,y2bar,xybar))
}

scoreMCE_theta = function(theta,ux,uy,N,Delta,s_dp,ep,n){
  synth = sdp_vec(ux,uy,N,Delta,theta,ep,n)
  synth = rbind(synth,s_dp)
  
  D_synth = depth.Mahalanobis(synth,synth)  
  r = rank(D_synth, ties.method="max")[R+1]
  s = r + D_synth[R+1]
  return(-s)
}

scoreMCE_nuisance = function(nuisance,beta1,ux,uy,N,Delta,s_dp,ep,n){
  theta = c(beta1,nuisance)  
  return(scoreMCE_theta(theta,ux,uy,N,Delta,s_dp,ep,n))
}


accept = function(beta1,ux,uy,N,Delta,s_dp,ep,n){
  rank_continuous = -scoreMCE_theta(c(beta1,nuisance),ux,uy,N,Delta,s_dp,ep,n)
  print(c(beta1, rank_continuous))
  if(rank_continuous >= (floor((alpha)*(R+1))+1))
    return(TRUE)
  
  opt = optim(par = nuisance, fn = scoreMCE_nuisance, lower=c(-5,-5,0.001,0.001), upper=c(5,5,5,5),
              method="L-BFGS-B",
              beta1=beta1,ux=ux,uy=uy,N=N,Delta=Delta,s_dp=s_dp,ep=ep,n=n
  )
  return((-opt$value)>=(floor((alpha)*(R+1))+1))
}





######################

# Figure 5 (left subfigure) is using ep=1 (here ep is the mu in mu-GDP)
# Appendix Figure 1 is using ep=0.5, 2, 0.2, 5, 0.1, 10, 0.05, 0.01
n_list <- c(100,200,300,400,500,1000,2000,5000)
beta1_list <- c(1, 0.8, 0.6, 0.4, 0.2, 0.1, 0)
ep_list <- c(1, 0.5, 2, 0.2, 5, 0.1, 10, 0.05, 0.01)
for (ep in ep_list) {
  
  nrow <- length(n_list)
  ncol <- length(beta1_list)
  power_table <- matrix(nrow = nrow, ncol = ncol)
  
  for (j in c(1:ncol)) {
    for (i in c(1:nrow)) {
      n <- n_list[i] 
      beta1_true <- beta1_list[j]
      
      beta[2] = beta1_true
      theta = c(beta1_true,beta[1],mu,tau,sa)
      nuisance = c(beta[1],mu,tau,sa)
      results <- foreach(
        rep = 1:reps,
        .combine = 'cbind',
        .packages=c('ddalpha'),
        .options.snow=opts
      ) %dopar% {
        # {rep=1
        set.seed(rep)
        ux = matrix(rnorm(1*n),ncol=n,nrow=1)
        uy = matrix(rnorm(1*n),ncol=n,nrow=1)
        N = matrix(rnorm(1*5),ncol=5,nrow=1)
        s_dp = sdp_vec(ux,uy,N,Delta,theta,ep,n)
        
        ### repros:
        
        
        ux = matrix(rnorm(R*n),ncol=n,nrow=R)
        uy = matrix(rnorm(R*n),ncol=n,nrow=R)
        N = matrix(rnorm(R*5),ncol=5,nrow=R)
        
        if(accept(beta1=0,ux,uy,N,Delta,s_dp,ep,n))
          rejected = 0
        else
          rejected = 1
        rejected
      }
      # print(rejected)
      power_table[i,j] <- mean(results)
      power_table2 <- data.frame(power_table)
      rownames(power_table2) <- n_list
      colnames(power_table2) <- beta1_list
      
      txt_name <- paste("./results/Repro_LR_HT/Repro_LR_HT", "-clamp_", Delta, 
                        "-gdp_ep=", ep, "-alpha=", alpha,
                        "-X_mu=", mu, "-tau=", tau,
                        "-beta=(", beta[1], ", b1)",
                        "-sa=", sa,
                        "-R=", R, "-reps=", reps, ".csv", sep='')
      print("")
      print(txt_name)
      print(t(power_table2))
      write.table(t(power_table2), txt_name, sep=",", row.names=TRUE, col.names=TRUE)
    }
  }
}



######################
# Figure 6 (first row of subfigures)
n_list <- c(100,200,300,400,500,1000,2000,5000)
ep <- 1
beta1_list <- c(1, 0)
Delta_list <- c(0.5, 0.8, 1.0, 1.5, 2, 5, 10)
for (Delta in Delta_list) {
  # for (ep in ep_list) {
  
  nrow <- length(n_list)
  ncol <- length(beta1_list)
  power_table <- matrix(nrow = nrow, ncol = ncol)
  
  for (j in c(1:ncol)) {
    for (i in c(1:nrow)) {
      n <- n_list[i]
      beta1_true <- beta1_list[j]
      
      beta[2] = beta1_true
      theta = c(beta1_true,beta[1],mu,tau,sa)
      nuisance = c(beta[1],mu,tau,sa)
      results <- foreach(
        rep = 1:reps,
        .combine = 'cbind',
        .packages=c('ddalpha'),
        .options.snow=opts
      ) %dopar% {
        set.seed(rep)
        ux = matrix(rnorm(1*n),ncol=n,nrow=1)
        uy = matrix(rnorm(1*n),ncol=n,nrow=1)
        N = matrix(rnorm(1*5),ncol=5,nrow=1)
        s_dp = sdp_vec(ux,uy,N,Delta,theta,ep,n)
        
        ### repros:
        
        
        ux = matrix(rnorm(R*n),ncol=n,nrow=R)
        uy = matrix(rnorm(R*n),ncol=n,nrow=R)
        N = matrix(rnorm(R*5),ncol=5,nrow=R)
        
        if(accept(beta1=0,ux,uy,N,Delta,s_dp,ep,n))
          rejected = 0
        else
          rejected = 1
        rejected
      }
      power_table[i,j] <- mean(results)
      power_table2 <- data.frame(power_table)
      rownames(power_table2) <- n_list
      colnames(power_table2) <- beta1_list
      
      txt_name <- paste("./results/compare_clamp/Repro_LR_HT", "-clamp_", Delta,
                        "-gdp_ep=", ep, "-alpha=", alpha,
                        "-X_mu=", mu, "-tau=", tau,
                        "-beta=(", beta[1], ", b1)",
                        "-sa=", sa,
                        "-R=", R, "-reps=", reps, ".csv", sep='')
      print("")
      print(txt_name)
      print(t(power_table2))
      write.table(t(power_table2), txt_name, sep=",", row.names=TRUE, col.names=TRUE)
    }
  }
}