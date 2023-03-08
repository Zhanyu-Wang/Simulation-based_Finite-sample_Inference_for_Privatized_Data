library(geometry) 
library(ddalpha)
library(tictoc)
library(rmutil)

################################################################
for(depth_type in c('mahalanobis', 'halfspace', 'simplicial', 'spatial')){
  n=100
  seed = rnorm(n)
  N = rlaplace(n=2,m=0,s=1)
  mu=1
  sa=1
  ep=1
  u=3
  l=0
  
  sdp_fun = function(seed,N,l,u,ep){
    n=length(seed)
    data = sa*seed+mu
    data_clamp = pmax(pmin(data,u),l)
    s1 = mean(data_clamp) + (u-l)/(n*ep)*N[1]
    s2 = var(data_clamp) + (u-l)^2/((n)*ep)*N[2]
    return(c(s1,s2))
  }
  s_dp = sdp_fun(seed,N,l,u,ep)
  
  
  
  sdp_vec = function(seed,N,l,u,ep,sa,mu){
    #seed should be Rxn
    #N should be Rx2
    n=dim(seed)[2]
    data = sa*seed+mu
    data_clamp = pmax(pmin(data,u),l)
    s1 = apply(data_clamp,1,mean) + (u-l)/(n*ep)*N[,1]
    s2 = apply(data_clamp,1,var) + (u-l)^2/((n)*ep)*N[,2]
    return(cbind(s1,s2))
  }
  
  score = function(theta,seed,N,n,s_dp,alpha){
    mu=theta[1]
    sa = theta[2]
    synth_sim = sdp_vec(seed,N,l,u,ep,sa=sa,mu=mu)
    synth = rbind(synth_sim,matrix(s_dp,ncol=2))
    set.seed(100)
    
    if(depth_type == 'mahalanobis'){
      D_synth = depth.Mahalanobis(synth,synth)
    }
    if(depth_type == 'halfspace'){
      D_synth = depth.halfspace(synth,synth)
    }
    if(depth_type == 'simplicial'){
      D_synth = depth.simplicial(synth,synth)
    }
    if(depth_type == 'spatial'){
      D_synth = depth.spatial(synth,synth)
    }
    
    rank = rank(D_synth)[R+1]
    return(-rank-D_synth[R+1])
  
  }
  accept = function(sa,mu,seed,N,n,s_dp,alpha){
    
    rank = -score(theta=c(mu,sa),seed,N,n,s_dp,alpha)
  
    if(rank>=(floor((alpha)*(R+1))+1)){
      return(TRUE)
    }
    else
      return(FALSE)
  
  }
  
  
  
  
  
  
  
  getBoundary = function(angle,theta_hat,s_dp,seed,N,n,ep,alpha){
    theta_0 = rep(0,2)
    theta_0[1]=theta_hat[1]### mu
    theta_0[2]=theta_hat[2]### sa
    
    theta_1 = rep(0,2)
    theta_1[1]=theta_hat[1] + 1*cos(angle*(2*pi))
    theta_1[2]=theta_hat[2] + 1*sin(angle*(2*pi))### what if negative??
    
    tol=10^-4
    while(sqrt(sum((theta_0-theta_1)^2))>tol){
      proposal = (theta_0+theta_1)/2
    
      if(accept(proposal[2],proposal[1],seed,N,n,s_dp,alpha))
        theta_0=proposal
      else{
        theta_1=proposal
      }
    }
    return(theta_1)
  }
  
  
  ### Plotting the confidence set
  n=100
  
  mu=1
  sa=1
  ep=1
  u=3
  l=0
  alpha=.05
  
  R = 200### for CI
  
  
  degrees = 500
  mu_boundary = rep(0,degrees)
  sa_boundary = rep(0,degrees)
  
  s_dp = c(1,.75)
  set.seed(235987)
  
  seed = matrix(rnorm(n*R,m=0,s=1),ncol=n,nrow=R)
  N = matrix(rnorm(R*2,m=0,s=1),ncol=2,nrow=R)
  tic()
  
  opt = optim(par=c(1,1),fn=score,lower=c(-5,0),upper=c(5,5),method="L-BFGS-B",seed=seed,N=N,n=n,s_dp=s_dp,alpha=alpha)
  theta_hat=opt$par
  
  pb <- txtProgressBar(max = degrees, style = 3)
  for(i in 1:degrees){
    setTxtProgressBar(pb, i)
    boundary = getBoundary(angle=i/degrees,theta_hat,s_dp,seed,N,n,ep,alpha)
    mu_boundary[i]=boundary[1]
    sa_boundary[i] = boundary[2]
  }
  toc()
  
  #### choose one depth
  pdf(paste("repro_normal_", depth_type, ".pdf", sep=""), width=5, height=5)  
  
  plot(mu_boundary,sa_boundary,type="l",
       xlab=expression(mu),ylab=expression(sigma),xlim=c(.4,1.22),ylim=c(.6,2.2))
  polygon(mu_boundary,sa_boundary,col="#1b98e0")
  dev.off() 

  
  area = sum((1/2)*((mu_boundary-theta_hat[1])^2+(sa_boundary-theta_hat[2])^2)*(2*pi/degrees))
  print(area)
}

