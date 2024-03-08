n = 100
ep = 1
l = 0
theta_true = 10
tol = 10^-4
reps = 1000
alpha = 0.05
R = 200

###################################################
# automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doSNOW",
  "ddalpha",
  "parallelly",
  "rmutil"
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
pb <- txtProgressBar(max = reps, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

start_time <- Sys.time()
#################################################

sdp_vec = function(e,N,theta,u,l,ep){
  data = qexp(e,rate=1/theta)
  data_clamp = pmax(pmin(data,u),l)
  S = apply((data_clamp),1,mean)
  return(S+((u-l)/(n*ep))*N)
}

scoreMCE = function(theta, e, N, u, l, ep, s_dp){
  synth_s = sdp_vec(e, N, theta, u, l, ep)
  synth = c(synth_s, s_dp)
  D_synth = 1 / (1 + (synth - mean(synth_s))^2/var(synth_s))
  r = rank(D_synth, ties.method="max")[R+1]
  s = r + D_synth[R+1]
  return(-s)
}


accept <- function(optim_par, s_dp, e, N, u, l, ep, search_left_bound, search_right_bound) {
  # check whether the proposal (optim_par) is acceptable
  proposed_result <- scoreMCE(optim_par, e, N, u, l, ep, s_dp)
  if ((-proposed_result) >= (floor((alpha)*(R+1))+1)){
    return(optim_par)
  }

  opt <-
    optim(
      par = optim_par,
      fn = scoreMCE,
      lower = c(search_left_bound),
      upper = c(search_right_bound),
      method = "L-BFGS-B",
      e = e,
      N = N,
      u = u,
      l = l,
      ep = ep,
      s_dp = s_dp
    )
  if ((-opt$value) >= (floor((alpha)*(R+1))+1)){
    return(opt$par)
  } else {
    # cannot find any param that is acceptable
    return(NA)
  }
}

getCI <- function(s_dp, e, N, u, l, ep) {
  theta_lower_bound <- 1e-4
  theta_upper_bound <- 1000
  # we try to find a theta_m in the confidence interval
  theta_m <- 10

  optim_par <- theta_m
  search_left_bound <- theta_lower_bound
  search_right_bound <- theta_upper_bound
  optim_par <- accept(optim_par, s_dp, e, N, u, l, ep, search_left_bound, search_right_bound)
  if (is.na(optim_par)){
    print("failed to find a point in CI")
    return(0)
  }
  theta_accepted_middle_val <- optim_par
  
  ## get CI left end
  left_lower <- theta_lower_bound
  left_upper <- theta_accepted_middle_val - tol
  optim_par_new <- optim_par
  while (left_upper - left_lower > tol){
    left_middle <- (left_lower + left_upper) / 2
    # check left half
    search_left_bound <- left_lower
    search_right_bound <- left_middle
    optim_par_new <- search_left_bound # (search_left_bound + search_right_bound) / 2
    optim_par_new <- accept(optim_par_new, s_dp, e, N, u, l, ep, search_left_bound, search_right_bound)
    if (is.na(optim_par_new)){
      # try it again with different init
      optim_par_new <- (search_left_bound + search_right_bound) / 2
      optim_par_new <- accept(optim_par_new, s_dp, e, N, u, l, ep, search_left_bound, search_right_bound)

      if (is.na(optim_par_new)){
        # exclude the left half
        left_lower <- search_right_bound
        optim_par_new <- optim_par
      }
      else {
        # find something in the left half
        left_upper <- optim_par_new - tol
      }
    }
    else{
      # find something in the left half
      left_upper <- optim_par_new - tol
    }
  }
  
  ## get CI right end
  right_lower <- theta_accepted_middle_val + tol
  right_upper <- theta_upper_bound
  optim_par_new <- optim_par
  while (right_upper - right_lower > tol) {
    right_middle <- (right_lower + right_upper) / 2
    # check right half
    search_left_bound <- right_middle
    search_right_bound <- right_upper
    optim_par_new <- search_right_bound # (search_left_bound + search_right_bound) / 2
    optim_par_new <- accept(optim_par_new, s_dp, e, N, u, l, ep, search_left_bound, search_right_bound)
    if (is.na(optim_par_new)){
      # try it again with different init
      optim_par_new <- optim_par
      optim_par_new <- (search_left_bound + search_right_bound) / 2
      optim_par_new <- accept(optim_par_new, s_dp, e, N, u, l, ep, search_left_bound, search_right_bound)

      if (is.na(optim_par_new)){
        # exclude the left half
        right_upper <- search_left_bound
        optim_par_new <- optim_par
      }
      else {
        # find something in the right half
        right_lower <- optim_par_new + tol
      }
    }
    else{
      # find something in the right half
      right_lower <- optim_par_new + tol
    }
  }
  
  return(c(left_lower, right_upper))
}


for(u in c(10, 20, 50, 100)){
# for(u in c(10)){
  results <- foreach(
    rep = 1:reps, 
    .combine = 'cbind',
    .packages=c('ddalpha', 'rmutil'),
    .options.snow=opts
  ) %dopar% {  
    # {rep=1
    set.seed(rep)
    
    e_true = matrix(runif(n*1,min=0,max=1),nrow=1,ncol=n)
    N_true = rlaplace(1,m=0,s=1)
    s_dp = sdp_vec(e_true, N_true, theta_true, u, l, ep)
    
    e = matrix(runif(n*R,min=0,max=1),nrow=R,ncol=n)
    N = rlaplace(R,m=0,s=1)
    
    CI <- getCI(s_dp, e, N, u, l, ep)
    current_covered <- (CI[1] <= theta_true & CI[2] >= theta_true)
    current_width <- CI[2] - CI[1]
    rbind(current_covered, current_width, CI[1], CI[2])
  }
  
  description_str <- paste("repro_exp-theta=", theta_true, 
                          "_alpha=", alpha,
                          "_ep=", ep,
                          "_R=", R,
                          "_n=", n,
                          "_clamp_l_u=", l,
                          "_", u,
                          "_reps=", reps,
                          sep='')
  txt_name <- paste(description_str,
                    ".results.csv", sep='')
  write.table(t(results), txt_name, sep=",", row.names=FALSE, col.names=TRUE)

  covered <- results[1,]
  width <- results[2,]
  result_stat <- data.frame("coverage"=mean(covered),
                            "se_coverage"=sqrt(mean(covered) * (1 - mean(covered)) / reps),
                            "mean_width"=mean(width),
                            "se_width"=sqrt(var(width) / reps),
                            "median_width"=median(width),
                            "Total time"=Sys.time() - start_time)

  txt_name <- paste(description_str,
                    ".stat.csv", sep='')
  write.table(t(result_stat), txt_name, sep=",", row.names=TRUE, col.names=TRUE)
  print(result_stat)
}
stopCluster(cl)
