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

# show progress bar
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)


start_time <- Sys.time()

clamp = function(x,a,b){
  return(pmin(pmax(x,a),b))
}


truncMean = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  return(mu+(dnorm(alpha)-dnorm(beta))*sa/Z)
}

truncVar = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  return((sa^2)*(1+(alpha*dnorm(alpha)-beta*dnorm(beta))/Z-((dnorm(alpha)-dnorm(beta))/Z)^2))
}

clampMean = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  #eye = Z==0
  #return((1-eye)*
  if(Z==0)
    return(a*pnorm(alpha)+b*(1-pnorm(beta)))
  else
    return(Z*truncMean(mu,sa,a,b)+a*pnorm(alpha)+b*(1-pnorm(beta)))
}

clampE2 = function(mu,sa,a,b){
  alpha = (a-mu)/sa
  beta = (b-mu)/sa
  Z = pnorm(beta)-pnorm(alpha)
  if(Z==0)
    EX2 = (a^2*pnorm(alpha)+b^2*(1-pnorm(beta)))
  else
    EX2 = (a^2*pnorm(alpha)+b^2*(1-pnorm(beta))+
             Z*(truncVar(mu,sa,a,b)+(truncMean(mu,sa,a,b))^2))
  return(EX2)
}

clampVar = function(mu,sa,a,b){
  EX2 = clampE2(mu,sa,a,b)
  return(EX2-(clampMean(mu,sa,a,b))^2)
}


clampXY = function(beta,mu,sa,tau,Delta){
  integrand = function(x){
    base = dnorm(x,mu,tau)
    # factor1 = pmax(pmin(x,Delta),-Delta)
    # factor2 = clampMean(beta[1]+beta[2]*x,sa,-Delta,Delta)
    # result = factor1*factor2*base
    factor12 = clampMean(x*(beta[1]+beta[2]*x),x^2*sa,-Delta^2,Delta^2)
    result = factor12*base
    if(is.nan(result)){
      print("Nan in ClampXY")
      return(0)
    }
    else
      return(result)
  }
  integrand_par = function(x){
    sapply(x,integrand)
  }
  #values = seq(-3,3,by=.01)
  #plot(values, integrand_par(values),type="l")
  return(integrate(integrand_par,lower=-Inf,upper=Inf)$value)
}

clampXY2 = function(beta,mu,sa,tau,Delta){
  integrand = function(x){
    base = dnorm(x,mu,tau)
    # factor1 = (pmax(pmin(x,Delta),-Delta))^2
    # factor2 = clampE2(beta[1]+beta[2]*x,sa,-Delta,Delta)
    # #factor2 = clampVar(beta[1]+beta[2]*x,sa,-Delta,Delta)+(clampMean(beta[1]+beta[2]*x,sa,-Delta,Delta))^2
    # result = factor1*factor2*base
    factor12 = clampE2(x*(beta[1]+beta[2]*x),x^2*sa,-Delta^2,Delta^2)
    result = factor12*base
    if(is.nan(result)){
      print("Nan in clampXY2")
      return(0)
    }
    else
      return(result)
  }
  integrand_par = function(x){
    sapply(x,integrand)
  }
  #values = seq(-3,3,by=.01)
  #plot(values, integrand_par(values),type="l")
  return(integrate(integrand_par,lower=-Inf,upper=Inf)$value)
}


clampXYvar = function(beta,mu,sa,tau,Delta){
  return(clampXY2(beta,mu,sa,tau,Delta)-(clampXY(beta,mu,sa,tau,Delta))^2)
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
  # xybar = apply((xc*yc),1,mean)+2*Delta^2/(n*ep)*N[,5]
  xybar = apply((xyc),1,mean)+2*Delta^2/(n*ep)*N[,5]
  
  return(data.frame(xbar,ybar,x2bar,y2bar,xybar))
}

obj_nuisance = function(nuisance,beta1,s_dp,Delta){
  beta0=nuisance[1]
  mu=nuisance[2]
  tau=nuisance[3]
  sa = nuisance[4]
  
  Ex = clampMean(mu,tau,-Delta,Delta)
  Ey = clampMean(beta0+beta1*mu,sqrt(beta1^2*tau^2+sa^2),-Delta,Delta)
  Ex2 = clampE2(mu,tau,-Delta,Delta)
  Ey2 = clampE2(beta0+beta1*mu,sqrt(beta1^2*tau^2+sa^2),-Delta,Delta)
  #Exy = clampXY(c(beta0,beta1),mu,sa,tau,Delta)
  return((s_dp[1]-Ex)^2+
           (s_dp[2]-Ey)^2+
           (s_dp[3]-Ex2)^2+
           (s_dp[4]-Ey2)^2)
}

nuisance_hat = function(beta1,s_dp,Delta,init=NULL){
  if(is.null(init))
    init=c(1,1,1,1)
  opt = optim(par=init,obj_nuisance,lower=c(-5,-5,0.001,0.001),upper=c(5,5,5,5),
              method="L-BFGS-B",
              beta1=beta1,s_dp=s_dp,Delta=Delta)
  return(opt$par)
}

getRank = function(theta,ux,uy,N,Delta,s_dp,ep,n){
  synth = sdp_vec(ux,uy,N,Delta,theta,ep,n)
  synth = rbind(synth,s_dp)
  
  D_synth = depth.Mahalanobis(synth,synth)
  r = rank(D_synth)[R+1]
  return(r)
}

scoreMCE_theta = function(theta,ux,uy,N,Delta,s_dp,ep,n){
  t = (floor((alpha)*(R+1))+1)
  synth = sdp_vec(ux,uy,N,Delta,theta,ep,n)
  synth = rbind(synth,s_dp)
  
  D_synth = depth.Mahalanobis(synth,synth)
  D_sort = sort(D_synth)
  # return(D_sort[t]-D_synth[R+1])
  r = sum(D_sort<=D_synth[R+1])
  s = r + D_synth[R+1]
  return(-s)
}

scoreMCE_nuisance = function(nuisance,beta1,ux,uy,N,Delta,s_dp,ep,n){
  theta = c(beta1,nuisance)
  
  return(scoreMCE_theta(theta,ux,uy,N,Delta,s_dp,ep,n))
}


accept = function(beta1,ux,uy,N,Delta,s_dp,ep,n){
  rank0 = getRank(c(beta1,nuisance),ux,uy,N,Delta,s_dp,ep,n)
  if(rank0>=(floor((alpha)*(R+1))+1))
    return(TRUE)
  
  n_hat = nuisance_hat(beta1=beta1,s_dp,Delta,init=nuisance)
  rank1 = getRank(c(beta1,n_hat),ux,uy,N,Delta,s_dp,ep,n)
  if(rank1>=(floor((alpha)*(R+1))+1))
    return(TRUE)
  
  ### optim fails with larger n...
  opt = optim(par =n_hat,fn = scoreMCE_nuisance,lower=c(-5,-5,0.001,0.001),upper=c(5,5,5,5)
              ,method="L-BFGS-B",
              beta1=beta1,ux=ux,uy=uy,N=N,Delta=Delta,s_dp=s_dp,ep=ep,n=n
  )
  return((-opt$value)>=(floor((alpha)*(R+1))+1))
}


getCI = function(ux,uy,N,Delta,s_dp,ep,n){
  
  if(accept(beta1=beta1_true,ux,uy,N,Delta,s_dp,ep,n))
    theta_m=beta1_true
  else{
    opt = optim(par=theta,fn=scoreMCE_theta,lower=c(-10,-10,-10,0.001,0.001),upper=c(10,10,10,10,10),
                method="L-BFGS-B",ux=ux,uy=uy,N=N,Delta=Delta,s_dp=s_dp,ep=ep,n=n)
    if(opt$value<=0)#getRank(opt$par,ux,uy,N,Delta,s_dp,ep,n)>=(floor((alpha)*(R+1))+1))
      theta_m=opt$par[1]
    else{
      print("Cannot find a point in CI")
      break
    }
  }
  ### binary search:
  theta_0=-10
  theta_1=10
  
  ##get lower bound:
  lower =theta_0
  middle=theta_m
  tol=10^-3
  while(middle-lower>tol){
    # print(c(lower,middle))
    proposal = (middle+lower)/2
    
    if(accept(proposal,ux,uy,N,Delta,s_dp,ep,n))
      middle=proposal
    else{
      lower=proposal
    }
  }
  #lower
  #accept(lower)
  #accept(middle)
  
  
  ##get upper bound:
  middle =theta_m
  upper=theta_1
  tol=10^-3
  while(upper-middle>tol){
    # print(c(middle,upper))
    proposal = (upper+middle)/2
    #synth = M(e,N,theta=proposal,U,ep)
    #interval=emp.hpd(synth)
    #interval = quantile(synth,c(.025,.975))
    if(accept(proposal,ux,uy,N,Delta,s_dp,ep,n))
      middle=proposal
    else{
      upper=proposal
    }
  }
  
  return(c(lower,upper))
}





######################

# Figure 4 (left subfigure) is using ep=1 (here ep is the mu in mu-GDP)
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
# Figure 5 (first row of subfigures)
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