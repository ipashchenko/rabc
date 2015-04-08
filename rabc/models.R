# Function that reads matrix of (bl, v0, tb) values and returns vectors of fluxes
flux <- function(x) {
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

flux_ell <- function(x) {
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

# Function must accept vector of model parameters and return summary statistics
# x - vector of parameters of population distribution (mu_logv0, std_logv0, mu_logtb, std_logtb)
# Example:
# x <- c(-0.43, 0.94, 11.8, 1.1)

model1 <- function(data, borders) {
  function(x) {
    fractions <- vector(length=length(borders)-1)
    for (i in seq(1, length(borders)-1)) {
      subdata<-subset(data, bl>borders[i] & bl<borders[i+1])
      # Size of sample in current bin of baselines
      size<-length(subdata$bl)
      # Generate sample of sources
      mu_logv0<-x[1]
      std_logv0<-x[2]
      mu_logtb<-x[3]
      std_logtb<-x[4]
      logv0 = rnorm(size, mu_logv0, std_logv0)
      logtb = rnorm(size, mu_logtb, std_logtb)
      fluxes<-flux(cbind(subdata$bl, 10**(logv0), exp(logtb)))
      # Count detections
      detections<-subset(subdata, fluxes>5.*subdata$s_thr)
      fractions[i]<-as.double(length(detections$bl))/size
    }
    return(fractions)
  }
}

# Function must accept vector of model parameters and return summary statistics
# x - vector of parameters of population distribution:
# (mu_logv0, std_logv0, mu_logtb, std_logtb, beta_e)
# Example:
# x <- c(-0.43, 0.94, 11.8, 2.1, 5.)
model2 <- function(data, borders) {
  function(x) {
    fractions <- vector(length=length(borders)-1)
    for (i in seq(1, length(borders)-1)) {
      subdata<-subset(data, bl>borders[i] & bl<borders[i+1])
      # Size of sample in current bin of baselines
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
      fluxes<-flux_ell(cbind(subdata$bl, exp(logv0), 10**(logtb), e, dfi))
      # Count detections
      detections<-subset(subdata, fluxes>5.*subdata$s_thr)
      fractions[i]<-as.double(length(detections$bl))/size
    }
    return(fractions)
  }
}

# Model with each single source having the same parameters
# even if it is observed on different baselines
# function must accept vector of model parameters and return summary statistics
# x - vector of parameters of population distribution:
# (mu_logv0, std_logv0, mu_logtb, std_logtb, beta_e)
# Example:
# x <- c(-0.43, 0.94, 11.8, 2.1, 5.)
model3 <- function(data, borders) {
  function(x) {
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
      # Size of sample in current bin of baselines
      n_obs_in_bin <-length(subdata$bl)
      dfi = runif(n_obs_in_bin, 0, pi/2)
      # Sources in current baseline bins
      sources_in_bin <- subdata$source
      v0 <- vector(length=n_obs_in_bin)
      tb <- vector(length=n_obs_in_bin)
      e <- vector(length=n_obs_in_bin)
      for (j in seq(1, n_obs_in_bin)) {
        v0[j] = exp(sources_params$logv0[sources_params$source == sources_in_bin[j]])
        tb[j] = 10**(sources_params$logtb[sources_params$source == sources_in_bin[j]])
        e[j] = sources_params$e[sources_params$source == sources_in_bin[j]]
      }
    
      fluxes<-flux_ell(cbind(subdata$bl, v0, tb, e, dfi))
      # Count detections
      detections<-subset(subdata, fluxes>5.*subdata$s_thr)
      fractions[i]<-as.double(length(detections$bl))/n_obs_in_bin
    }
    return(fractions)
  }
}