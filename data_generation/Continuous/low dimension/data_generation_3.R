library(MASS)

## Simulation setting ##

# path.cc <- "/Users/fanyiyang/Desktop/JAS\ revision/data/"
path.cc <- "/nas/longleaf/home/fyy2025/survival_tree_revision/lin_300_50/data/"
#path.cc <- "C:/Users/Haolin Li/Desktop/Dissertation/05_project2/02_simulation/05_RF/01_data/"
gamma.weibull = 0.7
p = 10
n = 300
nsim = 500
tau = 0.5
cmax = 20
cenc = rep(0, nsim)

## Data generation ##

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
#The autocorrelation structure of an AR(1) process describes the relationship between 
#the current observation and its previous 1 observation

# #tree structure
# true.lin.pred <- 10*(0.3*(X[,1]>0) -0.2*(X[,2]<0)+0.3*(X[,4]<0)*(X[,6]<0)-0.4*(X[,7]>0)*(X[,3]>0))
# #linear covariate
# true.lin.pred <- 10*(0.3*X[,1] -0.4*X[,2]+0.3*X[,4] -0.1*X[,6]-0.4*X[,7]+0.25*X[,3]+0.15*X[,8]+0.3*X[,9]-0.2*X[,10])
# #nonlinear covariate
# true.lin.pred <- 10*(0.3*X[,1] -0.2*X[,2]*X[,3]+0.3*X[,4]*X[,6]-0.4*X[,7]^3+0.1*X[,3]*X[,5]-0.25*X[,8]-0.3*X[,9]*X[,10]*X[,1])
# 10*(0.3*X[,11] -0.4*X[,12]+0.3*X[,14] -0.1*X[,16]-0.4*X[,17]+0.25*X[,13]+0.15*X[,18]+0.3*X[,19]-0.2*X[,20])

# low-dimensional tree-based DGP
for (j in 1:nsim){
  cat(j)
  X = as.matrix(mvrnorm(n = n, matrix(0, nrow = p, ncol = 1), ar1_cor(p, 0.6))) #mvrnorm: multivariate normal dist
  #mu: a vector of 0 with length p-- means of variables
  #Sigma: a positive-definite symmetric matrix specifying the covariance matrix of the variables.
  true.lin.pred <- 10*(0.3*X[,1] -0.4*X[,2]+0.3*X[,4] -0.1*X[,6]-0.4*X[,7]+0.25*X[,3]+0.15*X[,8]+0.3*X[,9]-0.2*X[,10])
  v <- runif(n=n)
  Tlat <- (- log(v) / (1* exp(true.lin.pred)))^(1 / gamma.weibull)
  # C = Tlat
  # censor_index = sample(n, n/2, replace = FALSE)
  # censor_time = runif(n/2, min = 0, max = cmax) # censoring times
  # for (i in 1:length(censor_index)){
  #   C[censor_index[i]]=censor_time[i]
  # }
  C = runif(n, min = 0, max = cmax) # censoring times
  end <- rep(tau, n)
  time <- pmin(Tlat, C, end) # follow-up times
  status <- as.numeric(Tlat == time) #event indicators
  dat = data.frame(X, status, time)
  write.csv(dat,file=paste0(path.cc, 'train', j, '.csv'), row.names = F)
  cenc[j] = table(dat$status)[1]/n
}

mean(cenc) # output censoring proportion

for (j in 1:nsim){
  cat(j)
  X = as.matrix(mvrnorm(n = n, matrix(0, nrow = p, ncol = 1), ar1_cor(p, 0.6)))
  true.lin.pred <- 10*(0.3*X[,1] -0.4*X[,2]+0.3*X[,4] -0.1*X[,6]-0.4*X[,7]+0.25*X[,3]+0.15*X[,8]+0.3*X[,9]-0.2*X[,10])
  v <- runif(n=n)
  Tlat <- (- log(v) / (1* exp(true.lin.pred)))^(1 / gamma.weibull)
  # C = Tlat
  # censor_index = sample(n, n/2, replace = FALSE)
  # censor_time = runif(n/2, min = 0, max = cmax) # censoring times
  # for (i in 1:censor_index){
  #   C[censor[i]]=censor_time[i]
  # }
  C = runif(n, min = 0, max = cmax) # censoring times
  end <- rep(tau, n)
  time <- pmin(Tlat, C, end) # follow-up times
  status <- as.numeric(Tlat == time) #event indicators
  test = data.frame(X, status, time)
  write.csv(test,file=paste0(path.cc, 'test', j, '.csv'), row.names = F)
  cenc[j] = table(test$status)[1]/n
}
mean(cenc)

