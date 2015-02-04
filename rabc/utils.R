library(gdata)
col.names <- c('source', 'bl', 's_thr', 'status')
setwd("/home/ilya/code/rabc/rabc")
data = read.table('data_to_R.txt', col.names=col.names)

borders <- c(2., 5., 10., 15., 20., 30.)

det_fractions_in_bsl_ranges <- function(data, borders){
  fractions <- vector(length=length(borders)-1)
  for (i in seq(1, length(borders)-1)) {
    statuses<-subset(data, bl>borders[i] & bl<borders[i+1], status)
    y<-subset(statuses, statuses == 'y')
    print(y)
    fractions[i] <- as.double(length(y$status))/length(statuses$status)
  }
  return(fractions)
}

#def flux_(b, v0, tb):
# Function that reads matrix of (bl, v0, tb) values and returns vectors of fluxes
flux<-function(x){
  # Transform baselines from ED to m
  x[ , 1] = x[ , 1] * 12742. * 10. ** 3
  # Transform fluxes from Jy to W / (m**2 * Hz)
  x[ , 2] = x[ , 2] * 10 ** (-26)
  k = 1.38 * 10 ** (-23)
  pi = 3.14159
  result = x[ , 2] * exp(-pi * x[ , 1] ** 2. * x[ , 2] / (2. * k * x[ , 3]))
  # Back to Jy
  result = result * 10 ** 26
  return(result)
}

#def flux(b, v0, tb, e, dfi)
flux_ell<-function(x){
  # Transform baselines from ED to m
  x[ , 1] = x[ , 1] * 12742. * 10. ** 3
  # Transform fluxes from Jy to W / (m**2 * Hz)
  x[ , 2] = x[ , 2] * 10 ** (-26)
  k = 1.38 * 10 ** (-23)
  pi = 3.14159
  result = x[ , 2] * exp(-pi * x[ , 1] ** 2. * x[ , 2] * (x[ , 4] * cos(x[ , 5]) ** 2. + x[ , 4] ** (-1) * sin(x[ , 5]) ** 2.) / (2. * k * x[ , 3]))
  # Back to Jy
  result = result * 10 ** 26
  return(result)
}

# function must accept vector of model parameters and return summary statistics
# x - vector of parameters of population distribution (mu_logv0, std_logv0, mu_logtb, std_logtb)
# example
# x<-c(-0.43, 0.94, 28.8, 2.1)
model1<-function(x){
  fractions <- vector(length=length(borders)-1)
  for (i in seq(1, length(borders)-1)) {
    subdata<-subset(data, bl>borders[i] & bl<borders[i+1])
    # Size of sample in current baseline bin
    size<-length(subdata$bl)
    # Generate sample of sources
    mu_logv0<-x[1]
    std_logv0<-x[2]
    mu_logtb<-x[3]
    std_logtb<-x[4]
    logv0 = rnorm(size, mu_logv0, std_logv0)
    logtb = rnorm(size, mu_logtb, std_logtb)
    fluxes<-flux(cbind(subdata$bl, exp(logv0), exp(logtb)))
    # Count detections
    detections<-subset(subdata, fluxes>5.*subdata$s_thr)
    fractions[i]<-as.double(length(detections$bl))/size
  }
  return(fractions)
}


# function must accept vector of model parameters and return summary statistics
# x - vector of parameters of population distribution:
# (mu_logv0, std_logv0, mu_logtb, std_logtb, beta_e)
# example
# x<-c(-0.43, 0.94, 28.8, 2.1, 5.)
model2<-function(x){
  fractions <- vector(length=length(borders)-1)
  for (i in seq(1, length(borders)-1)) {
    subdata<-subset(data, bl>borders[i] & bl<borders[i+1])
    # Size of sample in current baseline bin
    size<-length(subdata$bl)
    # Generate sample of sources
    mu_logv0<-x[1]
    std_logv0<-x[2]
    mu_logtb<-x[3]
    std_logtb<-x[4]
    beta_e<-x[5]
    logv0 = rnorm(size, mu_logv0, std_logv0)
    logtb = rnorm(size, mu_logtb, std_logtb)
    e = rbeta(size, 5, beta_e)
    dfi = runif(size, 0, pi/2)
    fluxes<-flux_ell(cbind(subdata$bl, exp(logv0), exp(logtb), e, dfi))
    # Count detections
    detections<-subset(subdata, fluxes>5.*subdata$s_thr)
    fractions[i]<-as.double(length(detections$bl))/size
  }
  return(fractions)
}

# as model2 but for multicore
model2m<-function(x){
  set.seed(x[1])
  fractions <- vector(length=length(borders)-1)
  for (i in seq(1, length(borders)-1)) {
    subdata<-subset(data, bl>borders[i] & bl<borders[i+1])
    # Size of sample in current baseline bin
    size<-length(subdata$bl)
    # Generate sample of sources
    mu_logv0<-x[2]
    std_logv0<-x[3]
    mu_logtb<-x[4]
    std_logtb<-x[5]
    beta_e<-x[6]
    logv0 = rnorm(size, mu_logv0, std_logv0)
    logtb = rnorm(size, mu_logtb, std_logtb)
    e = rbeta(size, 5, beta_e)
    dfi = runif(size, 0, pi/2)
    fluxes<-flux_ell(cbind(subdata$bl, exp(logv0), exp(logtb), e, dfi))
    # Count detections
    detections<-subset(subdata, fluxes>5.*subdata$s_thr)
    fractions[i]<-as.double(length(detections$bl))/size
  }
  return(fractions)
}

# TODO: Implement model with each single source having the same parameters
# even if it is observed on different baselines
# function must accept vector of model parameters and return summary statistics
# x - vector of parameters of population distribution:
# (mu_logv0, std_logv0, mu_logtb, std_logtb, beta_e)
# example
# x<-c(-0.43, 0.94, 28.8, 2.1, 5.)
model3<-function(x){
  fractions <- vector(length=length(borders)-1)
  for (i in seq(1, length(borders)-1)) {
    subdata<-subset(data, bl>borders[i] & bl<borders[i+1])
    # Size of sample in current baseline bin
    size<-length(subdata$bl)
    # Generate sample of sources
    mu_logv0<-x[1]
    std_logv0<-x[2]
    mu_logtb<-x[3]
    std_logtb<-x[4]
    beta_e<-x[5]
    logv0 = rnorm(size, mu_logv0, std_logv0)
    logtb = rnorm(size, mu_logtb, std_logtb)
    e = rbeta(size, 5, beta_e)
    dfi = runif(size, 0, pi/2)
    fluxes<-flux_ell(cbind(subdata$bl, exp(logv0), exp(logtb), e, dfi))
    # Count detections
    detections<-subset(subdata, fluxes>5.*subdata$s_thr)
    fractions[i]<-as.double(length(detections$bl))/size
  }
  return(fractions)
}


###################################################################
####################EasyABC########################################
###################################################################
library("EasyABC")
library("abc")
# Define priors
prior1=list(c("normal", -0.43, 0.15), c("normal", 0.94, 0.15), c("unif", 26, 32), c("unif", 0, 7.5))
prior2=list(c("normal", -0.43, 0.15), c("normal", 0.94, 0.15), c("unif", 26, 32), c("unif", 0, 7.5), c("unif", 1, 30))
# Data summary statistics
sum_stat_obs=c(0.68085106, 0.43043478, 0.20297030, 0.09770115, 0.02116402)
set.seed(1)

# Doing rejection sampling
ABC_rej<-ABC_rejection(model=model1, prior=prior1, nb_simul=100000,
                       summary_stat_target=sum_stat_obs, tol=0.005)

abc_out<-abc(sum_stat_obs, ABC_rej$param, ABC_rej$stats, tol=0.5,
             method="loclinear")

# Plot results
hist(abc_out$adj.values[, 1], plot="True")
hist(ABC_rej$param[, 3], plot="True", main="Histogram of log(Tb)", xlab="log(Tb)", ylab="P(log(Tb) | data)", freq=FALSE)
hist(abc_out$adj.values[, 3], plot="True", main="Histogram of log(Tb)", xlab="log(Tb)", ylab="P(log(Tb) | data)", freq=FALSE)
hist(ABC_rej$param[, 5], plot="True", main="Histogram of e", xlab="e", ylab="P(e | data)", freq=FALSE)
hist(abc_out$adj.values[, 5], plot="True", main="Histogram of e", xlab="e", ylab="P(e | data)", freq=FALSE)

# Doing PMC sampling
# This setup (20&0.02) results in 400s of computations
ABC_seq<-ABC_sequential(method="Drovandi", model=model2, prior=prior2, nb_simul=100,
                        summary_stat_target=sum_stat_obs, tolerance_tab=0.01)
ABC_seq<-ABC_sequential(method="Drovandi", model=model2m, prior=prior2, nb_simul=50,
                        summary_stat_target=sum_stat_obs, tolerance_tab=0.05, n_cluster=4, use_seed=TRUE)
abc_seq_out<-abc(sum_stat_obs, ABC_seq$param, ABC_seq$stats, tol=0.5,
                 method="loclinear")

###################################################################
###########Simulate survey#########################################
###################################################################

# First, give "TRUE" parameters of population
mu_logv0 = -0.8
std_logv0 = 0.3
mu_logtb = 29.0
std_logtb = 0.7
alpha_e = 33
beta_e = 50.

# Create ``n_s`` sources from population
n_s = length(unique(data$source))
sources = levels(data$sources)
set.seed(1)
logv0 = rnorm(n_s, mu_logv0, std_logv0)
logtb = rnorm(n_s, mu_logtb, std_logtb)
e = rbeta(n_s, alpha_e, beta_e)

# Create container for new data
new_data<-data
# Filling new_data with params
new_data$angles <- vector(length=length(new_data$bl))
new_data$logv0 <- vector(length=length(new_data$bl))
new_data$logtb <- vector(length=length(new_data$bl))
new_data$e <- vector(length=length(new_data$bl))
s=0
for (i in seq(along=sources)) {
  indxs<-which(new_data$source==new_data$source[i])
  s<-s+length(indxs)
  
  for (j in seq(along=indxs)) {
    new_data$e[indxs[j]]<-e[i]
    new_data$logv0[indxs[j]]<-logv0[i]
    new_data$logtb[indxs[j]]<-logtb[i]
    new_data$angles[indxs[j]]<-runif(1, 0, pi/2)
  }
}




fluxes<-flux_ell(cbind(subdata$bl, exp(logv0), exp(logtb), e, dfi))
# Count detections
detections<-subset(subdata, fluxes>5.*subdata$s_thr)
fractions[i]<-as.double(length(detections$bl))/size