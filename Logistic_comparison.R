library(ggplot2)
library(tidyverse)
####################################################################

reps_ori <- 1000

R <- 200
n = 500
alpha=.05
ep=1
shape1=1/2
shape2=1/2
beta0=1/2
beta1=beta1_true=2

n_list = c(100, 200, 500, 1000, 2000)
ep_list = c(0.1, 0.3, 1, 3, 10)

df_tab_coverage = data.frame(matrix(rep(0, length(n_list) * length(ep_list)), nrow=length(n_list), ncol=length(ep_list)))
colnames(df_tab_coverage) = ep_list
rownames(df_tab_coverage) = n_list
df_tab_width = df_tab_coverage
df_tab_count = df_tab_coverage

coverage_list = c()
width_list = c()
for(n in n_list){
  for(ep in ep_list){
    filename2 = paste("./results/summary/ERM_Logistic",
                      "-shape1_", shape1, "-shape2_", shape2,
                      "-beta0_", beta0, "-beta_1", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps_ori, ".csv", sep="")
    if (!file.exists(filename2)) {
      print(filename2)
      next
    }
    CIs = read.csv(filename2)
    reps = min(1000, length(CIs$l))
    covered = rep(0, reps)
    width = rep(0, reps)
    time_list = rep(0, reps)
    for (rep in c(1:reps)){
      if(CIs$l[rep]<=beta1_true & CIs$h[rep]>=beta1_true)
        covered[rep] = TRUE
      width[rep] = CIs$h[rep] - CIs$l[rep]
    }

    df_tab_coverage[as.character(n), as.character(ep)] = mean(covered[1:rep])
    df_tab_count[as.character(n), as.character(ep)] = reps
    df_tab_width[as.character(n), as.character(ep)] = mean(width[1:rep])
    coverage_list = rbind(coverage_list, c(ep, as.character(n), mean(covered[1:rep]), sd(covered[1:rep])/sqrt(reps)))
    width_list = rbind(width_list, c(ep, n, mean(width[1:rep]), sd(width[1:rep])))
    # width_list = rbind(width_list, c(ep, n, mean(width[1:rep]), quantile(width[1:rep], c(0.25,0.75))))
  }
}

df_tab_coverage
write.csv(t(df_tab_coverage), "./results/ERM_logistic.csv")
df_tab_width
df_tab_count
df_coverage_list = data.frame(coverage_list)
ep_mapping = c('solid','dashed','dotted','longdash')
colnames(df_coverage_list) = c('ep', 'n', 'coverage', 'coverage_std')
df_width_list = data.frame(width_list)
colnames(df_width_list) = c('ep', 'n', 'width', 'width_std')
# colnames(df_width_list) = c('ep', 'n', 'width', 'width_quantile1', 'width_quantile3')
df_width_list$ep = factor(as.character(df_width_list$ep), levels = c("0.1", "0.3", "1", "3", "10"))
df_width_list$lt = df_width_list$ep


####################################################################

df_tab_coverage = data.frame(matrix(rep(0, length(n_list) * length(ep_list)), nrow=length(n_list), ncol=length(ep_list)))
colnames(df_tab_coverage) = ep_list
rownames(df_tab_coverage) = n_list
df_tab_width = df_tab_coverage
df_tab_count = df_tab_coverage

coverage_list = c()
width_list = c()
for(n in n_list){
  for(ep in ep_list){
    filename2 = paste("./results/summary/Repro_Logistic",
                      "-shape1_", shape1, "-shape2_", shape2,
                      "-beta0_", beta0, "-beta_1", beta1, "-R_", R, "-n_", n,
                      "-alpha_", alpha, "-ep_", ep, "-reps_", reps_ori, ".csv", sep="")
    # print(filename2)
    if (!file.exists(filename2)) {
      next
    }
    CIs = read.csv(filename2)
    # reps = length(CIs$V1)
    reps = min(1000, length(CIs$V1))
    covered = rep(0, reps)
    width = rep(0, reps)
    time_list = rep(0, reps)
    for (rep in c(1:reps)){
      if(CIs$V1[rep]<=beta1_true & CIs$V2[rep]>=beta1_true)
        covered[rep] = TRUE
      width[rep] = CIs$V2[rep] - CIs$V1[rep]
    }
    
    df_tab_coverage[as.character(n), as.character(ep)] = mean(covered[1:rep])
    df_tab_count[as.character(n), as.character(ep)] = reps
    df_tab_width[as.character(n), as.character(ep)] = mean(width[1:rep])
    coverage_list = rbind(coverage_list, c(ep, as.character(n), mean(covered[1:rep]), sd(covered[1:rep])/sqrt(reps)))
    width_list = rbind(width_list, c(ep, n, mean(width[1:rep]), sd(width[1:rep])))
  }
}

df_tab_coverage
write.csv(t(df_tab_coverage), "./results/Repro_logistic.csv")
df_tab_width
df_tab_count
df_coverage_list = data.frame(coverage_list)
ep_mapping = c('solid','dashed','dotted','longdash')
colnames(df_coverage_list) = c('ep', 'n', 'coverage', 'coverage_std')
df_width_list2 = data.frame(width_list)
colnames(df_width_list2) = c('ep', 'n', 'width', 'width_std')
df_width_list2$ep = factor(as.character(df_width_list2$ep), levels = c("0.1", "0.3", "1", "3", "10"))
df_width_list2$lt = df_width_list2$ep

##################################################


df_width_list$method = "DP-CI-ERM"
df_width_list2$method = "Repro Sample"
df_width_list3 <- rbind(df_width_list2, df_width_list)

override.linetype <- (c(1,2,3,4,5))
ggplot(data=df_width_list3, aes(x = n, y = width, group=ep, lt=lt, colour=factor(ep))) +
  geom_line(aes(linetype=lt)) +
  geom_errorbar(aes(ymin=pmax(0.31, width-width_std), ymax=width+width_std), width=.01,
                position=position_dodge(.01)) +
  scale_y_log10(limits = c(0.3, 26), breaks = c(0.1, 0.3,0.5, 1, 2,5,10, 20)) +
  scale_x_log10(limits = c(90, 2100), breaks = c(100,200,500, 1000, 2000)) +
    labs(title=NULL,
       x ="sample size", y = "CI width") +
  scale_linetype_manual(values=c('solid','dashed','dotted','dotdash','longdash'))+
  theme_bw() + theme(plot.title = element_text(hjust = 0), legend.position="top") + 
  guides(color=guide_legend(reverse = FALSE, title=expression(epsilon*"  "), nrow=1, 
                            override.aes = list(linetype = override.linetype)), 
         linetype="none", alpha="none") +
  facet_wrap(~ fct_rev(method), nrow=1)
ggsave("logistic_width_comparison.pdf", width=5, height=3)
