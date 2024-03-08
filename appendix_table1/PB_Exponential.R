n=100
ep = 1
u= 10
for(u in c(10, 20, 50, 100)){
  l=0
  theta_true=10
  set.seed(42)
  
  
  M_exp = function(e,N,theta,u,l,ep){
    data = qexp(e,rate=1/theta)
    data_clamp = pmax(pmin(data,u),l)
    S = apply((data_clamp),1,mean)
    return(S+((u-l)/(n*ep))*N)
  }
  
  
  alpha=.05
  reps = 1000
  cover_pivot = cover_db = rep(0,reps)
  width_pivot = width_db = rep(0,reps)
    
  pb <- txtProgressBar(max = reps, style = 3)
  for(i in 1:reps){
    setTxtProgressBar(pb, i)
    e_true = matrix(runif(n*1,min=0,max=1),nrow=1,ncol=n)
    N_true = rlaplace(1,m=0,s=1)
  
    s_dp = M_exp(e_true,N_true,theta_true,u,l,ep)
    s_dp
    
    R=1000
    e = matrix(runif(n*R,min=0,max=1),nrow=R,ncol=n)
    N = rlaplace(R,m=0,s=1)
    
    synth = M_exp(e,N,theta=max(s_dp,0),u,l,ep)
    qs = quantile(synth-s_dp,c(.025,.975))
    
    CI = c(max(s_dp-qs[2],0),s_dp-qs[1])
    if(CI[1]<=theta_true & theta_true<=CI[2])
      cover_pivot[i]=TRUE
    
  
    CI_db = quantile(synth-2*(mean(synth)-s_dp),c(.025,.975))
    if(CI_db[1]<=theta_true & theta_true<=CI_db[2])
      cover_db[i]=TRUE
    
    width_pivot[i] = CI[2]-CI[1]
    width_db[i] = CI_db[2]-CI_db[1]
  }
    
  print(c(u, mean(cover_pivot), sqrt(mean(cover_pivot)*(1-mean(cover_pivot))/reps)))
  print(c(u, mean(width_pivot), sqrt(var(width_pivot)/reps)))
  # mean(cover_db)
  # mean(width_db)
}
  