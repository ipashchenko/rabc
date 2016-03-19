set.seed(1)
band='c'
data_lookup = c("c"="data_to_R_Cband.txt", "l"="data_to_R_Lband.txt",
                "k"="data_to_R_Kband.txt")

load_data <- function(fname, path="/home/ilya/code/rabc/rabc") {
  library(gdata)
  col.names <- c('source', 'bl', 's_thr', 'status')
  setwd(path)
  data = read.table(fname, col.names=col.names)
  return(data)
}

create_tb <- function(prob, mus, sds, n) {
  components <- sample(1:2,prob=prob,size=n,replace=TRUE)
  tb_samples <- rnorm(n=n,mean=mus[components],sd=sds[components])
  return(tb_samples)
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

create_survey<-function(data, mu_logv0=-0.4, std_logv0=0.3, mu_logtb=12.0, std_logtb=0.7, alpha_e=33, beta_e=50,
                        prob=NaN) {
  # Set seed for reproducibility
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  set.seed(1)
  # Number of observations in survey
  n_obs = length(data$bl)
  # Create ``n_s`` sources from population with given parameters
  n_s = length(unique(data$source))
  sources = levels(data$source)
  logv0 = rnorm(n_s, mu_logv0, std_logv0)
  if (!is.nan(prob)) {
    logtb <- create_tb(prob, mu_logtb, std_logtb, n=n_s)
  }
  else {
    logtb = rnorm(n_s, mu_logtb, std_logtb)
  }
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
  new_data$status[new_data$fluxes > 5*new_data$s_thr] = 1.
  new_data$status[new_data$fluxes < 5*new_data$s_thr] = 0.
  # Delete some columns
  new_data$angles <- NULL
  new_data$logv0 <- NULL
  new_data$logtb <- NULL
  new_data$e <- NULL
  new_data$fluxes <- NULL

  return(new_data)
}

get.probs <- function(data) {
  my.logit <- glm(status ~ bl, data=data, family="binomial")
  probs <- 1./(1. + exp(-my.logit$coefficients[1]) - my.logit$coefficients[2] * data$bl)
  return(probs)
}

get.probs.gam <- function(data) {
  my.logit <- gam(status ~ lo(bl), data=data, family="binomial")
  return(my.logit$fitted.values)
}

model22 <- function(data) {
  probs <- get.probs(data)
  function(x) {
    new_data <- create_survey(data, mu_logv0=x[1], std_logv0=x[2],
                              mu_logtb=x[3], std_logtb=x[4], alpha_e=5.,
                              beta_e=x[5])
    new_probs <- get.probs(new_data)
    probs_diff <- new_probs - probs
    sum.stat <- sqrt(sum(probs_diff**2))
    # print(c("sum.stat for", x, " is ", sum.stat))
    return(sum.stat)
  }
}

data = load_data(data_lookup[band], path='/home/ilya/code/rabc/rabc/')
fake_data = create_survey(data, mu_logtb = c(12, 13.5), std_logtb = c(0.3, 0.3),
                          prob = c(0.85, 0.15))

sum_stat_obs <- 0.
library("EasyABC")
library("abc")
prior=list(c("normal", -0.43, 0.15),
         c("normal", 0.94, 0.15),
         c("unif", 7., 16.),
         c("unif", 0, 3.5),
         c("unif", 1, 30))
model <- model22(fake_data)
ABC_rej<-ABC_rejection(model=model, prior=prior, nb_simul=50000,
                       summary_stat_target=sum_stat_obs, tol=0.025,
                       progress_bar=TRUE)
# abc_out<-abc(sum_stat_obs, ABC_rej$param, ABC_rej$stats, tol=0.5,
#              method="neuralnet", numnet=50, sizenet=2)

