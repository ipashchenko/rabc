library(gdata)
col.names <- c('source', 'bl', 's_thr', 'status')
setwd("/home/ilya/code/rabc/rabc")
data = read.table('data_to_R.txt', col.names=col.names)

# Model that given parameters vector ``x`` returns vector of summary statistic
model <- function(x){
  
}

borders <- c(2., 5., 10., 15., 20.)
#test commit windows
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
  x[ , 1] = x[ , 1] * 12742. * 10. ** 3
  k = 1.38 * 10 ** (-23)
  pi = 3.14159
  result = x[ , 2] * exp(-pi * x[ , 1] ** 2. * x[ , 2] * 10 ** (-26) / (2. * k * x[ , 3]))
  return(result)
}

# example
x<-c(-0.43, 0.94, 28.8, 2.1)
# function must accept vector of model parameters and return summary statistics
# x - vector of parameters of population distribution (mu_logv0, std_logv0, mu_logtb, std_logtb)
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


###################################################################
####################EasyABC########################################
###################################################################

# Define priors
prior1=list(c("unif", -3, 1), c("unif", 0, 5), c("unif", 26, 32), c("unif", 0, 7.5))
# Data summary statistics
sum_stat_obs=c(0.68085106, 0.43043478, 0.20297030, 0.09770115)
tolerance=0.5
n=100
set.seed(1)
ABC_rej<-ABC_rejection(model=model1, prior=prior1, nb_simul=n,
                       summary_stat_target=sum_stat_obs, tol=tolerance)