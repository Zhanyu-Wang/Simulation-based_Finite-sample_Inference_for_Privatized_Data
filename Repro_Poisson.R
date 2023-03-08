library(tictoc)
width = rep(0,50)
lower_vec = rep(0,50)
upper_vec = rep(0,50)

n=100
ep = 1
l=0
theta_true=10

set.seed(12345)
e_true = matrix(runif(n*1,min=0,max=1),nrow=1,ncol=n)
N_true = rnorm(1,m=0,s=1)

R=1000
e = matrix(runif(n*R,min=0,max=1),nrow=R,ncol=n)
N = rnorm(R,m=0,s=1)


alpha=.05

pb <- txtProgressBar(max = 50, style = 3)
for(u in 1:50){
  setTxtProgressBar(pb, u)
  
  M = function(e,N,theta,u,l,ep){
    data = qpois(e,lambda=theta)
    data_clamp = pmax(pmin(data,u),l)
    S = apply((data_clamp),1,mean)
    return(S+((u-l)/(n*ep))*N)
  }
  s_dp = M(e_true,N_true,theta_true,u,l,ep)
  s_dp
  
  theta_m = theta_true
  synth = c(M(e,N,theta=theta_m,u,l,ep),s_dp)
  
  rank = rank(synth)[R+1]
  if(rank>=(floor((alpha)/2*(R+1))+1) & rank<=ceiling((1-alpha/2)*(R+1))){
    found_middle=TRUE
  }
  else {
    print("problem")
    break
  }
  
  ### binary search:
  theta_0=0
  theta_m=theta_true
  theta_1=1000
  
  ##get lower bound:
  lower =theta_0
  middle=theta_m
  tol=10^-4
  while(middle-lower>tol){
    proposal = (middle+lower)/2
    synth = c(M(e,N,theta=proposal,u,l,ep),s_dp)
    rank = rank(synth)[R+1]
    if(rank>=(floor((alpha)/2*(R+1))+1) & rank<=ceiling((1-alpha/2)*(R+1)))
      middle=proposal
    else{
      lower=proposal
    }
  }
  #lower
  
  ##get upper bound:
  middle =theta_m
  upper=theta_1
  tol=10^-4
  while(upper-middle>tol){
      proposal = (upper+middle)/2
      synth = c(M(e,N,theta=proposal,u,l,ep),s_dp)
      rank = rank(synth)[R+1]
      if(rank>=(floor((alpha)/2*(R+1))+1) & rank<=ceiling((1-alpha/2)*(R+1)))
        middle=proposal
      else{
        upper=proposal
      }
  }
  # print(upper)
  
  width[u] = upper-lower
  lower_vec[u] = lower
  upper_vec[u] = upper
}


pdf("widthPoisson.pdf",width=5,height=5)
plot(seq(1,50),(width),type="l",xlab="U",ylab="width",ylim=c(1,4))
points(seq(1,50),(width),pch=16)
dev.off()
### minimum is 1.372637 at U=14
which(width==min(width))
min(width)

pdf("CIPoisson.pdf",width=5,height=5)
lower_vec2 = c(0,lower_vec)
upper_vec2 = c(1000,upper_vec)
plot(c(0,0),c(lower_vec2[1],upper_vec2[1]),ylim=c(0,15),xlim=c(0,50),type="l",xlab="U",ylab="confidence intervals")
points(c(0,0),c(lower_vec2[1],upper_vec2[1]),pch="-")
for(i in 2:50){
  lines(c(i-1,i-1),c(lower_vec2[i],upper_vec2[i]))
  points(c(i-1,i-1),c(lower_vec2[i],upper_vec2[i]),pch="-")
}
dev.off()

#########################################################
###   coverage sanity check
#########################################################
set.seed(2134135)
n=100
ep = 1
u= 10
l=0
theta_true=10

alpha=.05

reps = 1000
coverage = 0

pb <- txtProgressBar(max = reps, style = 3)
for(i in 1:reps){
  setTxtProgressBar(pb, i)
  
  e_true = matrix(runif(n*1,min=0,max=1),nrow=1,ncol=n)
  N_true = rnorm(1,m=0,s=1)
  
  R=1000
  e = matrix(runif(n*R,min=0,max=1),nrow=R,ncol=n)
  N = rnorm(R,m=0,s=1)
  
  M = function(e,N,theta,u,l,ep){
    data = qpois(e,lambda=theta)
    data_clamp = pmax(pmin(data,u),l)
    S = apply((data_clamp),1,mean)
    return(S+((u-l)/(n*ep))*N)
  }
  s_dp = M(e_true,N_true,theta_true,u,l,ep)
  s_dp
  synth = c(M(e,N,theta=theta_true,u,l,ep),s_dp)
  
  rank = rank(synth)[R+1]
  if(rank>=(floor((alpha)/2*(R+1))+1) & rank<=ceiling((1-alpha/2)*(R+1)))
    coverage = coverage + 1/reps
}
coverage
sqrt(coverage*(1-coverage)/reps)

