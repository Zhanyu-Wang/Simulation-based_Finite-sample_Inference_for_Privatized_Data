###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha",
  "parallelly",
  "tictoc",
  "geometry"
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
reps = 1000### for coverage simulation
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
###################################################

sdp_fun = function(seed,N,ep,theta,n){
  #n=length(seed)
  B = qbinom(seed,size=n,prob=theta)
  s1 = (B) + (1/ep)*N[1]
  s2 = (n-(B))+ (1/ep)*N[2]
  return(c(s1,s2))
}
#s_dp = sdp_fun(seed,N,l,u,ep)



sdp_vec = function(seed,N,ep,theta,n){
  #seed should be Rxn
  #N should be Rx2
  #n=dim(seed)[2]
  #if(n>=1)
    B = qbinom(seed,size=n,prob=theta)
  #else 
  #  B=rep(0,R)
  s1 = B + (1/ep)*N[,1]
  s2 = (n-B)+ (1/ep)*N[,2]
  return(cbind(s1,s2))
}


getRank = function(n,theta,seed,N,s_dp,ep){
  synth = sdp_vec(seed,N,ep,theta=theta,n=n)
  synth = rbind(synth,s_dp)
  
  #depth.contours(synth, depth = "halfspace", main = "halfspace")
  #points(s_dp[1],s_dp[2])
  
  #D_synth = depth.sample(synth,synth)
  #D_synth = D_synth$Simpl
  
  D_synth = depth.Mahalanobis(synth,synth)
  ### could append s_dp at the end of 2nd synth
  #F_D = ecdf(D_synth)
  #D_sdp = depth.simplicial(data=synth,x=matrix(s_dp,ncol=2))
  #D_sdp = depth.Mahalanobis(data=synth,x=matrix(s_dp,ncol=2))
  rank = rank(D_synth)[R+1]
  
  #s = F_D(D_sdp)### needs to be bigger than 1-alpha (.05)
  return(rank)
}

test = function(s,theta,ep){
  n_hat = as.numeric(pmax(s[,1]+s[,2],1))
  #t = (as.numeric(s[,1])-n_hat*theta)/sqrt(n_hat*theta*(1-theta)+1/ep^2)
  #t = (as.numeric(s[,1])-n_hat*theta)/sqrt((theta^2+(1-theta)^2)*(n_hat*theta*(1-theta)+1/ep^2))
  t = (as.numeric(s[,1])-n_hat*theta)/sqrt((n_hat*theta*(1-theta)+(theta^2+(1-theta)^2)/ep^2))
  return(t)
}

scoreMCE= function(n,theta,seed,N,s_dp,ep,pivot=FALSE){
  n=floor(n)
  t = (floor((alpha)*(R+1))+1)
  synth = sdp_vec(seed,N,ep=ep,theta=theta,n=n)
  
  synth = rbind(synth,s_dp)
 
  if(pivot){
  tests = test(s=synth,theta,ep)
  #synth = matrix(synth[,1],nrow=R+1)
  synth = matrix(tests,nrow=(R+1))
  }
  D_synth = depth.Mahalanobis(synth,synth)
  D_sort = sort(D_synth)
  return(D_sort[t]-D_synth[R+1])#(rank-(R+1)/2)^2)
}



accept = function(theta,seed,N,s_dp,ep,pivot){
  ### optim fails with larger n...
  n_low = max(0,floor(sum(s_dp)-qnorm(.99995,0,sqrt(2)/ep)))
  n_high = max(0,ceiling(sum(s_dp)+qnorm(.99995,0,sqrt(2)/ep)))
  n_cand = seq(n_low,n_high)
  min_val = 1
  for(n in n_cand){
    val = scoreMCE(n,theta,seed,N,s_dp,ep,pivot)
    min_val = min(min_val,val)
  }
  #opt = optim(par = 50,fn = scoreMCE,lower=0,upper=1000,method="Brent",
  #            theta=theta,seed=seed,N=N,s_dp=s_dp,ep=ep)
  #opt = optim(par = 2,fn = score,lower=0,method="L-BFGS-B",mu=m)
  #rank = getRank(sa=opt$par,mu=mu,seed,N,n,s_dp,ep)
  #print(opt$value)
  #print(opt$par)
  if(min_val<=0)
    return(TRUE)
  else
    return(FALSE)
}


getCI = function(s_dp,seed,N,ep,pivot=FALSE){
  
  if(accept(theta=theta_true,seed,N,s_dp,ep,pivot))
    theta_m=theta_true
  else{
    #opt = optim(par=c(.5,100),fn=scoreMCE2,lower=c(0,0),upper=c(1,Inf),
    #            method="L-BFGS-B",seed=seed,N=N,s_dp=s_dp,ep=ep)
    theta_hat = min(max(s_dp[1],0)/max(sum(s_dp),1),1)
    #print(theta_hat)
    if(accept(theta=theta_hat,seed,N,s_dp,ep,pivot))
      theta_m=theta_hat
    
    else{
      print("failed to find a point in CI")
      return(c(0,0))
      #break
    }
    #opt2 = optim(par=1.5,fn=score,lower=-5,upper=10,method="Brent",seed=seed,N=N,n=n,s_dp=s_dp)
    #theta_m=opt2$par
  }
  ### binary search:
  theta_0=0
  theta_1=1
  
  ##get lower bound:
  lower =theta_0
  middle=theta_m
  tol=10^-4
  while(middle-lower>tol){
    #print(c(lower,middle))
    proposal = (middle+lower)/2
  
    if(accept(proposal,seed,N,s_dp,ep,pivot))
      middle=proposal
    else{
      lower=proposal
    }
  }
  #lower
  #accept(lower)
  #accept(middle)
  
  
  ##get upper bound:
  middle =lower
  upper=theta_1
  tol=10^-4
  while(upper-middle>tol){
    #print(c(middle,upper))
    proposal = (upper+middle)/2
    #synth = M(e,N,theta=proposal,U,ep)
    #interval=emp.hpd(synth)
    #interval = quantile(synth,c(.025,.975))
    if(accept(proposal,seed,N,s_dp,ep,pivot))
      middle=proposal
    else{
      upper=proposal
    }
  }
  
  return(c(lower,upper))
}

#getCI(s_dp,seed,N,n)



#lower = .64
#accept(.64,seed,N,n,s_dp)
#sigma = 1.448


### mini simulation
n_true=100

theta_true=.2
ep=1

alpha=.05

R = 200### for CI


covered = coveredP=rep(0,reps)
width = widthP=rep(0,reps)
tic()

results <- foreach(
  rep = 1:reps, 
  .combine = 'cbind',
  .packages=c('ddalpha'),
  .options.snow=opts
) %dopar% {
  set.seed(rep)
  setTxtProgressBar(pb, rep)
  ### get obs value
  seed = runif(1,0,1)
  N = rnorm(n=2,m=0,s=1)
  s_dp = sdp_fun(seed,N,ep,theta_true,n_true)
  
  seed = runif(R,0,1)
  N = matrix(rnorm(R*2,m=0,s=1),ncol=2,nrow=R)
  CI = getCI(s_dp,seed,N,ep)
  CI
  if(CI[2]==0)
    print(s_dp)
  if(CI[1]<=theta_true & CI[2]>=theta_true)
    covered[rep] = TRUE
  width[rep] = CI[2]-CI[1]
  #print(mean(covered[1:rep]))
  #print(mean(width[1:rep]))
  CIpivot = getCI(s_dp,seed,N,ep,pivot=TRUE)
  CIpivot
  if(CIpivot[2]==0)
    print(s_dp)
  if(CIpivot[1]<=theta_true & CIpivot[2]>=theta_true)
    coveredP[rep] = TRUE
  widthP[rep] = CIpivot[2]-CIpivot[1]
  print(mean(coveredP[1:rep]))
  print(mean(widthP[1:rep]))
  c(covered[rep], width[rep], coveredP[rep], widthP[rep])
}
stopCluster(cl)
toc()
covered <- results[1,]
width <- results[2,]
coveredP <- results[3,]
widthP <- results[4,]

print('Mahalanobis Depth coverage')
print(c(mean(covered), sqrt(mean(covered)*(1-mean(covered))/reps)))
print('Mahalanobis Depth width')
print(c(mean(width), sqrt(var(width)/reps)))

print('Approximate Pivot coverage')
print(c(mean(coveredP), sqrt(mean(coveredP)*(1-mean(coveredP))/reps)))
print('Approximate Pivot width')
print(c(mean(widthP), sqrt(var(widthP)/reps)))

