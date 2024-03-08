#!/usr/bin/env Rscript

reps <- 1000

R <- 200
n = 500
alpha=.05
ep=2
shape1=1/2
shape2=1/2
beta0=1/2
beta1=beta1_true=2

set.seed(42)
dir.create(file.path('.', 'results/logistic'), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path('.', 'results/summary'), showWarnings = FALSE, recursive = TRUE)

###
n_list = c(100, 200, 500, 1000, 2000)
ep_list = c(0.1, 0.3, 1, 3, 10)

pivot = TRUE

for(ep in ep_list){
  for(n in n_list){
    filename2 = paste("./results/summary/Repro_Logistic",
                      "-shape1_", shape1, "-shape2_", shape2,
                      "-beta0=", beta0, "-beta1=", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps, ".csv", sep="")
    print(filename2)
    
    CIs = c()
    for (rep in c(1:reps)){
      filename = paste("./results/logistic/Repro_Logistic",
                      "-shape1_", shape1, "-shape2_", shape2,
                      "-beta0=", beta0, "-beta1=", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps, "-", rep, ".csv", sep="")
      tryCatch({
        csv_file = read.csv(filename, header=TRUE)
        CI = csv_file$X1
        CIs = rbind(CIs, CI)
      },
      error=function(cond){print(filename)})
    }
    write.csv(CIs, filename2)
  }
}


