### parametric bootstrap
library(rmutil)
library(tictoc)
#library(truncnorm)
ep=1
n=100
m_true=1
s_true = 1
a = 0
b=3
clamp = function(x){
  return(pmin(pmax(x,a),b))
}
X_func = function(u,m,s){
  return(s*u+m)
}
R = 200
alpha=.05

coverage = coverage2=coverageS=coverageS2=coveragePB=coverageSPB=0
reps = 1000
width =width2= widthS=widthS2=widthPB=widthSPB=rep(0,reps)
set.seed(1000)

pb <- txtProgressBar(max = reps, style = 3)
for(r in 1:reps){
  setTxtProgressBar(pb, r)
  X_0 = rnorm(n,m_true,s_true)
  
  x_clamp = clamp(X_0)
  z_m = mean(x_clamp) + (rnorm(n=1,m=0,s=(b-a)/ep))/n
  z_s = var((x_clamp)) + (rnorm(n=1,m=0,s=(b-a)^2/ep)/(n)) #+ 4*(b-a)^2/(ep*(n-1))
  #z_s = var((x_clamp)) + rtruncnorm(n=1,mean=0,sd=(b-a)^2/(ep*(n-1)),a=-var((x_clamp))) #+ 4*(b-a)^2/(ep*(n-1))
  
  z_s_sqrt=sqrt(max(z_s,0))
  
  m_np = mean(X_0)
  s_np = sqrt(var(X_0))
  
  u = matrix(rnorm(R*n,m=0,s=1),nrow=R,ncol=n)
  w1 = rnorm(n=R,m=0,s=(b-a)/ep)
  w2 = rnorm(n=R,m=0,s=(b-a)^2/ep)
  m_star =m_star2= m_starPB=rep(0,R)
  s_star =s_star2= s_starPB=rep(0,R)
  

  for(i in 1:R){
    
    data = X_func(u[i,],m=z_m,s=sqrt(max(z_s,0)))
    data_clamp = clamp(data)
    z_star1 = mean(data_clamp)+(w1[i])/n
    z_star2 = var(data_clamp)+(w2[i])/(n)#+ 4*(b-a)^2/(ep*(n-1))
    #z_star2 = var(data_clamp)+rtruncnorm(n=1,mean=0,sd=(b-a)^2/(ep*(n-1)),a=-var((data_clamp)))

    m_star[i] = z_star1
    s_star[i] = sqrt(max(z_star2,0))
    m_star2[i]= 2*z_m-z_star1
    #s_star2[i] = 2*z_s-z_star2#max(2*z_s-z_star2,0)
    s_star2[i] = 2*z_s_sqrt-s_star[i]
    #data_np = X_func(u[i,],m=m_np,s=s_np)
    #m_starPB[i] = mean(data_np)
    #s_starPB[i] = var(data_np)
    
  }
  #toc()
  interval = quantile(m_star,c(alpha/2,1-alpha/2))
  interval2 = quantile(m_star2,c(alpha/2,1-alpha/2))
  #interval2 = c(z_m-qm[2],z_m-qm[1])
  #interval = quantile(m_star,c(.025,.975),na.rm=TRUE)
  if(interval[1]<=m_true & interval[2]>=m_true){
    coverage = coverage + 1/reps
  }
  width[r] = interval[2]-interval[1]
  if(interval2[1]<=m_true & interval2[2]>=m_true){
    coverage2 = coverage2 + 1/reps
  }
  width2[r] = interval2[2]-interval2[1]
  
  intervalS = quantile(s_star,c(alpha/2,1-alpha/2))
  #intervalS2=2*z_s-intervalS
  intervalS2 = quantile(s_star2,c(alpha/2,1-alpha/2))
  #intervalS2 = c(z_s-qs[2],z_s-qs[1])
  
  #intervalS = sqrt(pmax(intervalS,0))
  #intervalS2= sqrt(pmax(intervalS2,0))
  
  if(intervalS[1]<=s_true & intervalS[2]>=s_true){
    coverageS = coverageS + 1/reps
  }
  widthS[r] = intervalS[2]-intervalS[1]

  if(intervalS2[1]<=s_true & intervalS2[2]>=s_true){
    coverageS2 = coverageS2 + 1/reps
  }
  widthS2[r] = intervalS2[2]-intervalS2[1]
}
close(pb)

### Parametric Bootstrap (percentile) for mu (mean)
c(coverage, sqrt(coverage*(1-coverage)/reps))
c(mean(width), sqrt(var(width)/reps))

### Parametric Bootstrap (percentile) for sigma (standard deviation)
c(coverageS, sqrt(coverageS*(1-coverageS)/reps))
c(mean(widthS), sqrt(var(widthS)/reps))


### Parametric Bootstrap (simplified t) for mu (mean)
c(coverage2, sqrt(coverage2*(1-coverage2)/reps))
c(mean(width2), sqrt(var(width2)/reps))

### Parametric Bootstrap (simplified t) for sigma (standard deviation)
c(coverageS2, sqrt(coverageS2*(1-coverageS2)/reps))
c(mean(widthS2), sqrt(var(widthS2)/reps))
