# simulation settings
value_r = 100
upper_clamp <- 3
lower_clamp <- 0
n <- 100 # sample size
R_synthetic <- 200 # number of synthetic samples
ep <- 1 # privacy guarantee parameter
alpha <- .05 # confidence level
tol <- 10^-8

population_mu <- 1 # population mean
population_sigma <- 1 # population sd



################################################################
### Plotting the confidence set

dp_statistic = c(1, 0.75)  
d1_lower = -10
d1_upper = 10
d2_lower = tol
d2_upper = 10

set.seed(1)
data_randomness <- matrix(rnorm(n * R_synthetic), ncol = n, nrow = R_synthetic)
privacy_noises <- matrix(rnorm(R_synthetic * 2), ncol = 2, nrow = R_synthetic)

################################################################
for(depth_type in c('mahalanobis', 'halfspace', 'simplicial', 'spatial')){
  boundary_valid = read.csv(file=paste(depth_type, value_r, ep, "boundary.csv", sep='_'))
  boundary_valid = boundary_valid[,2:5]

  area = sum((boundary_valid[,3] - boundary_valid[,1]) * (boundary_valid[,4] - boundary_valid[,2]))
  print(c(depth_type, area))
  # ep=1, 100
  # 0.354115559299602 mahalanobis 0.5514926 1.1880036 0.7841644 1.5345536
  # 0.605506006334619 halfspace   0.312500  2.125000  0.750000  1.457031
  # 0.630322806823634 simplicial  0.2265625 2.1250000 0.7500000 1.4570313
  # 0.362314512180075 spatial     0.5568847 1.1849995 0.7767796 1.5360642

  # ep=0.5, 100
  # 0.766327478059884 mahalanobis 0.3612175  1.2197266  0.5502777  1.9074708
  # 0.953061589678131 halfspace   0.48370359 1.14062501 0.00000001 3.25000001
  # 0.964212564381321 simplicial  0.47957156 1.14062501 0.00000001 3.25000001
  # 0.776741354683335 spatial     0.3676376  1.2290956  0.5339965  1.9078484

  pdf(paste(depth_type, value_r, ep, "region.pdf", sep="_"), width = 5, height = 5)  
  plot(c(0.2, 2.2), c(0.7, 1.6), type = "n", xlab=expression(mu),ylab=expression(sigma)) # , main = area
  rect(boundary_valid[,1], 
       boundary_valid[,2], 
       boundary_valid[,3], 
       boundary_valid[,4], 
       border = NA,
       col = "#1b98e0")

  invisible(dev.off())
}
# stopCluster(cl)

################################################################



