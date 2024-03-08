library(TeachingDemos)
library(rmutil)
library(tictoc)
for(U in c(10, 20, 50, 100)){
  set.seed(42)
  
  M_exp = function(e,N,theta,u,l,ep){
    data = e*theta#qexp(e,rate=1/theta)
    data_clamp = pmax(pmin(data,u),l)
    S = apply((data_clamp),1,mean)
    return(S+((u-l)/(n*ep))*N)
  }
  
  inversion = function(q,theta,U,n,ep){
    integrand = function(t){
      i = 0+1i
      inv = exp(-i*t*q)
      main = ((1-exp(U*(i*(t/n)-1/theta)))/(1-i*(t/n)*theta)+
                (exp(i*(t/n)*U-U/theta)))^n
      Lap = 1/(1+(U/(n*ep))^2*t^2)
      ans = Im(inv*main*Lap)/t
      return(ans)
    }
    cdf = 1/2 - 1/pi*integrate(integrand,lower=0,upper=Inf)$value
    return(cdf)
  }
  
  objCDF = function(theta){
    return((inversion(s_dp,theta,U,n,ep)-1/2)^2)
  }
  
  inversionPDF = function(q,theta,U,n,ep){
    integrand = function(t){
      i = 0+1i
      inv = exp(-i*t*q)
      main = ((1-exp(U*(i*(t/n)-1/theta)))/(1-i*(t/n)*theta)+
                (exp(i*(t/n)*U-U/theta)))^n
      Lap = 1/(1+(U/(n*ep))^2*t^2)
      ans = (inv*main*Lap)
      return(ans)
    }
    integrandRe = function(t){
      return(Re(integrand(t)))
    }
    #integrandIm = function(t){
    #  return(Im(integrand(t)))
    #}
    pdfReal = integrate(integrandRe,lower=-Inf,upper=Inf)$value
    #pdfIm = integrate(integrandRe,lower=-Inf,upper=Inf)$value
    return(pdfReal)
  }
  
  
  
  objPDF = function(theta){
    return(-(inversionPDF(s_dp,theta,U,n,ep)))
  }
  
  getCI = function(s_dp,U,n,ep){
    if(inversion(s_dp,10,U,n,ep)>=.025 & inversion(s_dp,10,U,n,ep)<=.975)
      middle=10
    else{
      opt2 = optim(par=10,fn=objCDF,lower=0,upper=100,method="Brent")
      #print(opt2)
      middle=opt2$par
    }
    lower =0
    #middle=10
    tol=10^-6
    while(middle-lower>tol){
      #print(c(lower,middle))
      proposal = (middle+lower)/2
      #synth = M(e,N,theta=proposal,U,ep)
      #interval=emp.hpd(synth)
      #interval = quantile(synth,c(.025,.975))
      cdf = inversion(q=s_dp,theta=proposal,U,n,ep)
      if(cdf<=.975)
        middle=proposal
      else{
        lower=proposal
      }
    }
    
    
    #middle=10
    middle=lower
    upper = 100
    tol=10^-6
    while(upper-middle>tol){
      #print(c(middle,upper))
      proposal = (middle+upper)/2
      #synth = M(e,N,theta=proposal,U,ep)
      #interval=emp.hpd(synth)
      #interval = quantile(synth,c(.025,.975))
      cdf = inversion(q=s_dp,theta=proposal,U,n,ep)
      if(cdf>=.025)
        middle=proposal
      else{
        upper=proposal
      }
    }
    #lower
    #upper
    return(c(lower,upper))
  }
  
  
  #######################################
  ### CONFIDENCE INTERVALS
  ###
  
  n=100
  ep = 1
  theta_true=10
  
  R = 1000
  coverage = 0
  width = rep(0,R)
  
  pb <- txtProgressBar(max = R, style = 3)
  for(i in 1:R){
    setTxtProgressBar(pb, i)
    e_true = matrix(rexp(n*1,rate=1),nrow=1,ncol=n)
    N_true = rlaplace(1,m=0,s=1)
    s_dp = M_exp(e_true,N_true,theta_true,U,0,ep)
  
    CI = getCI(s_dp,U,n,ep)
    if(theta_true>=CI[1] & theta_true<=CI[2])
      coverage = coverage + 1/R
    #else{
    #  print(s_dp)
    #  print(CI)
    #}
    if(CI[2]>99)
      print(CI)
    width[i] = CI[2]-CI[1]
  }
  
  print(c(U, coverage, sqrt(coverage*(1-coverage)/R)))
  print(c(U, mean(width), sqrt(var(width)/R)))
}
