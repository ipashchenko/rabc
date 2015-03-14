library(gdata)
col.names <- c('source', 'bl', 's_thr', 'status')
setwd("/home/ilya/code/rabc/rabc")
data = read.table('data_to_R.txt', col.names=col.names)
# Read & Write objects example
# On odin
# saveRDS(ABC_rej3, 'ABC_rej3.rds')
# on calculon
# ABC_rej3 <- readRDS("ABC_rej3.rds")

# TODO: calculate fractions for different number of bins
fractions_in_bins <- function(some_data, n_bins) {
  borders <- vector(length=n_bins+1)
  d = (max(some_data$bl) - min(some_data$bl)) / n_bins
  for (i in seq(1, n_bins+1)) {
    borders[i] <- min(some_data$bl)+(i-1)*d
  }
  return(borders)
}


det_fractions_in_bsl_ranges <- function(some_data, borders){
  fractions <- vector(length=length(borders)-1)
  for (i in seq(1, length(borders)-1)) {
    statuses<-subset(some_data, bl>borders[i] & bl<borders[i+1], status)
    y<-subset(statuses, statuses == 'y')
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

# as model2 but for multiple cores
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
  # set seed for the same P(x) if x is the same
  fractions <- vector(length=length(borders)-1)
  # Create ``n_s`` sources from population
  n_s = length(unique(data$source))
  sources = levels(data$source)
  mu_logv0<-x[1]
  std_logv0<-x[2]
  mu_logtb<-x[3]
  std_logtb<-x[4]
  beta_e<-x[5]
  logv0 = rnorm(n_s, mu_logv0, std_logv0)
  logtb = rnorm(n_s, mu_logtb, std_logtb)
  e = rbeta(n_s, 5, beta_e)
  sources_params <- data.frame(source=sources, logv0=logv0, logtb=logtb, e=e)
  
  for (i in seq(1, length(borders)-1)) {
    # set seed for the same angles
    subdata<-subset(data, bl>borders[i] & bl<borders[i+1])
    # Size of sample in current baseline bin
    n_obs_in_bin <-length(subdata$bl)
    dfi = runif(n_obs_in_bin, 0, pi/2)
    # Sources in current baseline bins
    sources_in_bin <- subdata$source
    v0 <- vector(length=n_obs_in_bin)
    tb <- vector(length=n_obs_in_bin)
    e <- vector(length=n_obs_in_bin)
    for (j in seq(1, n_obs_in_bin)) {
      v0[j] = exp(sources_params$logv0[sources_params$source == sources_in_bin[j]])
      tb[j] = exp(sources_params$logtb[sources_params$source == sources_in_bin[j]])
      e[j] = sources_params$e[sources_params$source == sources_in_bin[j]]
    }
    
    fluxes<-flux_ell(cbind(subdata$bl, v0, tb, e, dfi))
    # Count detections
    detections<-subset(subdata, fluxes>5.*subdata$s_thr)
    fractions[i]<-as.double(length(detections$bl))/n_obs_in_bin
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
ABC_rej<-ABC_rejection(model=model1, prior=prior1, nb_simul=50000,
                       summary_stat_target=sum_stat_obs, tol=0.05,
                       progress_bar=TRUE)

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

myfunc <- function(x, i) {
  return((x[,i] - mean(x[,i])) / sd(x[,i]))
}

# TODO: Add standardisation based on variances of priors.
distances_k <- function(sample, k) {
  # Now msample[1,] - first element of sample
  msample <- matrix(sample, nrow=dim(sample)[1], ncol=dim(sample)[2])
  # Standardise 
  for (col in seq(1, ncol(msample))) {
    msample[,col]<-(msample[,col]-mean(msample[,col]))/sd(msample[,col])
  }
  distances <- vector(length = dim(sample)[1])
  for (i in seq(1, dim(sample)[1])) {
    element <- msample[i,]
    # Now melement[j,] - is i-th element for any j 
    melement <- t(matrix(rep(element, dim(msample)[1]), nrow=dim(sample)[2], ncol=dim(sample)[1]))
    diffs_sqrd <- apply(msample - melement, 2, function(x) x**2)
    sums_diffs_sqrd <- apply(diffs_sqrd, 1, sum)
    sqrt_sums_diffs_sqrd <- sqrt(sums_diffs_sqrd)
    # indx <- which(sqrt_sums_diffs_sqrd == sort(sqrt_sums_diffs_sqrd)[k])
    distances[i]<-sort(sqrt_sums_diffs_sqrd)[k] 
  }
  return(distances)
}

  
# Calculate entropy of sample
entropy <- function(sample, n_par, k=4) {
  # vector of k-th nearest neighbors for sample
  distances <- distances_k(sample, k)
  n <- dim(sample)[1]
  print(c("n", n))
  return(log(pi**(n_par/2)/gamma(n_par/2+1)) - digamma(k) + log(n) + (n_par/n)* sum(log(distances)))
}

# Select number of bins used for summary statistics
# for each n_bins_max=10 make n=5 rejection ABC samplings
# some_data here is data$bl
find_n_bins <- function(some_data, n_bins_max=7, n=5) {
  means_entropy_all = vector(length=n_bins_max)
  std_entropy_all = vector(length=n_bins_max)
  for (i in seq(1, n_bins_max)) {
    entropy_one = vector(length=n)
    # Claclulate data summary statistics
    borders <<- fractions_in_bins(some_data, i)
    sum_stat_observed <- det_fractions_in_bsl_ranges(some_data, borders)
    print("Using number of bins: ")
    print(i)
    print("Using borders: ")
    print(borders)
    print("Using summary statistics: ")
    print(sum_stat_observed)
    for (j in seq(1, n)) {
      print("Inside inner cycle, borders :")
      print(borders)
      print(c("Using", i, "number of bins", j, "time"))
      ABC_rej<-ABC_rejection(model=model2, prior=prior2, nb_simul=5000,
                             summary_stat_target=sum_stat_observed, tol=0.04,
                             progress_bar=TRUE)
      print("summary statistic from sample")
      print(ABC_rej$stats)
      print(length(sum_stat_observed))
      print(dim(ABC_rej$stats[2]))
      abc_out<-abc(sum_stat_observed, ABC_rej$param, ABC_rej$stats, tol=0.5, method="loclinear")
      entropy_one[j] = entropy_one[j] + entropy(abc_out$adj.values, 5, k=4)
    }
    means_entropy_all[i]=mean(entropy_one)
    std_entropy_all[i]=sd(entropy_one)
  }
  return(c(means_entropy_all, std_entropy_all))
}
 
###################################################################
###########Simulate survey#########################################
###################################################################

create_survey<-function(some_data, mu_logv0=-0.4, std_logv0=0.3, mu_logtb=29.0, std_logtb=0.7, alpha_e=33, beta_e=50) {

  set.seed(1)
  # Number of observations in survey
  n_obs = length(some_data$bl)
  # Create ``n_s`` sources from population
  n_s = length(unique(some_data$source))
  sources = levels(some_data$source)
  logv0 = rnorm(n_s, mu_logv0, std_logv0)
  logtb = rnorm(n_s, mu_logtb, std_logtb)
  e = rbeta(n_s, alpha_e, beta_e)

  # Create container for new data
  new_data<-some_data
  # Create new entries
  # TODO: Do i need this?
  new_data$angles <- vector(length=length(n_obs))
  new_data$logv0 <- vector(length=length(n_obs))
  new_data$logtb <- vector(length=length(n_obs))
  new_data$e <- vector(length=length(n_obs))
  new_data$fluxes <- vector(length=length(n_obs))
  # Filling new_data with params
  for (i in seq(along=sources)) {
    indxs<-which(new_data$source==sources[i])
    for (j in seq(along=indxs)) {
      new_data$e[indxs[j]]<-e[i]
      new_data$logv0[indxs[j]]<-logv0[i]
      new_data$logtb[indxs[j]]<-logtb[i]
      new_data$angles[indxs[j]]<-runif(1, 0, pi/2)
    }
  }

  # Calculate observed fluxes for sample
  new_data$fluxes<-flux_ell(cbind(new_data$bl, exp(new_data$logv0), exp(new_data$logtb), new_data$e, new_data$angles))
  # Update statuses
  new_data$status[new_data$fluxes > 5*new_data$s_thr] = 'y'
  new_data$status[new_data$fluxes < 5*new_data$s_thr] = 'n'

  return(new_data)
}


# Count detections
detections<-subset(subdata, fluxes>5.*subdata$s_thr)
fractions[i]<-as.double(length(detections$bl))/size

# Simulate from logit model
simulate.from.logr <- function(somedata, coefs) {
  require(faraway) # For accessible logit and inv-logit functions
  n = length(somedata$status)
  linear.part = coefs[1] + data$b * coefs[2] + data$d * coefs[3]
  probs = ilogit(linear.part)
  y = rbinom(n, size=1, prob=probs)
  return(y)
}

# Simulate from logistic fitted model and re-fit both logistic and GAM
delta.deviance.sim <- function(somedata, logistic.model) {
  y.new = simulate.from.logr(somedata, logistic.model$coefficients)
  GLM.dev = glm(y.new ~ b + d, data=somedata, family="binomial")$deviance
  GAM.dev = gam(y.new ~ lo(b) + lo(d), data=somedata, family="binomial")$deviance
  return(GLM.dev - GAM.dev)
}

#plot adj
plot(ABC_rej3$param[,3], sqrt(rowSums((ABC_rej3$stats )**2)), xlab="log(Tb)", ylab="S")
lf <- locfit(sqrt(rowSums((ABC_rej3$stats )**2))  ~  ABC_rej3$param[,3])
plot(lf, c(27,29.5), add=TRUE, lwd=5)
lvl <- sqrt(sum(sum_stat_obs**2))
lines(c(26, 30), c(lvl, lvl), lty=3)
lines(c(26, 30), c(lvl+0.1, lvl+0.1))
lines(c(26, 30), c(lvl-0.1, lvl-0.1))

# Save data to txt
write.table(ABC_rej3$param, "abc3_all.txt", sep="\t")
write.table(ABC_rej2$param, "abc2_all.txt", sep="\t")
write.table(ABC_rej1$param, "abc1_all.txt", sep="\t")
