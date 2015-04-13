# Function that given data frame with ``source``, ``bl``, ``s_thr``, ``status``
# columns return data frame with simulated survey.
create_survey<-function(data, mu_logv0=-0.4, std_logv0=0.3, mu_logtb=12.0, std_logtb=0.7, alpha_e=33, beta_e=50) {
  # Set seed for reproducibility
  set.seed(1)
  # Number of observations in survey
  n_obs = length(data$bl)
  # Create ``n_s`` sources from population with given parameters
  n_s = length(unique(data$source))
  sources = levels(data$source)
  logv0 = rnorm(n_s, mu_logv0, std_logv0)
  logtb = rnorm(n_s, mu_logtb, std_logtb)
  e = rbeta(n_s, alpha_e, beta_e)
  
  # Create container for new data
  new_data <- data
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
  new_data$fluxes<-flux_ell(cbind(new_data$bl, exp(new_data$logv0), 10**(new_data$logtb), new_data$e, new_data$angles))
  # Update statuses
  new_data$status[new_data$fluxes > 5*new_data$s_thr] = 'y'
  new_data$status[new_data$fluxes < 5*new_data$s_thr] = 'n'
  # Delete some columns
  new_data$angles <- NULL
  new_data$logv0 <- NULL
  new_data$logtb <- NULL
  new_data$e <- NULL
  new_data$fluxes <- NULL
  
  return(new_data)
}