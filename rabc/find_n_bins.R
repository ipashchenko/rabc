library("EasyABC")
library("abc")
source("/home/ilya/code/rabc/rabc/enthropy.R")
source("/home/ilya/code/rabc/rabc/data_load.R")
source("/home/ilya/code/rabc/rabc/models.R")
source("/home/ilya/code/rabc/rabc/simulate_survey.R")


get_prior_std <- function(prior_list) {
  prior_std <- vector(length=length(prior_list))
  for (i in seq(1, length(prior_list))) {
    prior <- prior_list[[i]]
    if (prior[1:1] == "normal") {
      std <- as.numeric(prior[3:3]) 
    }
    else if (prior[1:1] == "unif") {
      std <- (as.numeric(prior[3:3]) - as.numeric(prior[2:2])) / sqrt(12.)
    }
    prior_std[i] <- std
  }
  return(prior_std)
}

# Select number of bins used for summary statistics using enthropy
# for each number of bins from ``1`` to ``n_bins_max make`` make n=5 rejection
# ABC samplings
# data here is data$bl
find_n_bins_enthropy <- function(band, n_bins_max=7, n=10, path_to_data='/home/ilya/code/rabc/rabc/',
                                 prior=list(c("normal", -0.43, 0.15),
                                            c("normal", 0.94, 0.15),
                                            c("unif", 7., 16.),
                                            c("unif", 0, 3.5),
                                            c("unif", 1, 30)),
                                 nb_simul=50000, tol=0.05, standard="eaton") {
  # Load data
  data_lookup = c("c"="data_to_R_Cband.txt", "l"="data_to_R_Lband.txt",
                  "k"="data_to_R_Kband.txt")
  data = load_data(data_lookup[band], path=path_to_data)
  
  # Calculate std of priors for standardisation
  prior_std_list <- get_prior_std(prior)
  
  means_entropy_all = vector(length=n_bins_max)
  std_entropy_all = vector(length=n_bins_max)
  # Circle for number of bins in summary statistics
  for (i in seq(1, n_bins_max)) {
    entropy_one = vector(length=n)
    # Calculate data summary statistics
    bins <<- calc_bins(data, i)
    sum_stat_observed <- calc_fractions_in_bins(data, bins)
    print("Using number of bins: ")
    print(i)
    print("Using borders: ")
    print(bins)
    print("Using summary statistics: ")
    print(sum_stat_observed)
    for (j in seq(1, n)) {
      print("Inside inner cycle, borders :")
      print(bins)
      print(c("Using", i, "number of bins", j, "time"))
      cur_model <- model2(data, bins)
      ABC_rej<-ABC_rejection(model=cur_model, prior=prior, nb_simul=nb_simul,
                             summary_stat_target=sum_stat_observed, tol=tol,
                             progress_bar=TRUE)
      print("summary statistic from sample")
      print(ABC_rej$stats)
      print(length(sum_stat_observed))
      print(dim(ABC_rej$stats[2]))
      #abc_out<-abc(sum_stat_observed, ABC_rej$param, ABC_rej$stats, tol=0.5, method="loclinear")
      #entropy_one[j] = entropy_one[j] + entropy(abc_out$adj.values, 5, k=4)
      entropy_one[j] = entropy_one[j] + entropy(ABC_rej$param, 5, k=4, prior_std_list=prior_std_list)
    }
    means_entropy_all[i]=mean(entropy_one)
    std_entropy_all[i]=sd(entropy_one)
  }
  return(c(means_entropy_all, std_entropy_all))
}

# Function that first, using ME chosen S_ME (that is n_bins_ME) find n_obs = nb_simul * tol_ME / 2 parameter vectors
# and generate n_obs new data sets from each parameter vector. Then for each number of bins (1 to n_bins_max) for each
# data set it calculates RSSE and average them to MRSEE. Finally, it returns MRSEE for each number of bins
find_n_bins_MRSSE <- function(band, n_bins_ME, nb_simul_ME=10000, tol_ME=0.05,
                              n_bins_max=10, path_to_data='/home/ilya/code/rabc/rabc/',
                              prior=list(c("normal", -0.43, 0.15),
                                         c("normal", 0.94, 0.15),
                                         c("unif", 7., 16.),
                                         c("unif", 0, 3.5),
                                         c("unif", 1, 30)),
                              nb_simul_MRSSE=5000, tol_MRSSE=0.04,
                              standard=TRUE, use_correction=TRUE, n_obs=50) {
  # Vector to store MRSSE for each bin
  MRSSE_output <- vector(length=n_bins_max)
  
  # Load data
  data_lookup = c("c"="data_to_R_Cband.txt", "l"="data_to_R_Lband.txt",
                  "k"="data_to_R_Kband.txt")
  data = load_data(data_lookup[band], path=path_to_data)
  
  # Calculate std of priors for standardisation
  prior_std_list <- get_prior_std(prior)
  
  # Calculate bins for `best` ME number of bins
  bins <- calc_bins(data, n_bins_ME)
  # First, generate n_obs parameters
  cur_model <- model2(data, bins)
  sum_stat_observed <- calc_fractions_in_bins(data, bins)
  ABC_rej<-ABC_rejection(model=cur_model, prior=prior, nb_simul=nb_simul_ME,
                         summary_stat_target=sum_stat_observed, tol=tol_ME,
                         progress_bar=TRUE)
  abc_out<-abc(sum_stat_observed, ABC_rej$param, ABC_rej$stats, tol=0.5,
               method="neuralnet", numnet=50, sizenet=2)
  
  if (use_correction==FALSE) {
    sample = ABC_rej$param
  }
  else if (use_correction==TRUE) {
    sample = abc_out$adj.values
  }
  
  # Choose only n_obs number of parameter vectors from simulated data
  sample <- sample[sample(nrow(sample),size=n_obs,replace=FALSE),]
  
  # TODO: Use scale(sample, scale=as.vector(prior_std_list))
  # Standardize #############################
  #  if (is.null(standard)) {
  #    for (col in seq(1, ncol(sample))) {
  #      sample[,col]<-(sample[,col]-mean(sample[,col]))/sd(sample[,col])
  #    } 
  #  }
  # Keep it unstandardized to create new data sets
  params <- sample[]
  
  if (standard==TRUE) {
    # Standardize a-la Ewan et al. 
    # prior_std_list - list of std of priors for parameters
    # (0.15, 0.15, 6./sqrt(12.), 7.5/sqrt(12.), 29./sqrt(12.))
    for (col in seq(1, ncol(sample))) {
      sample[,col]<-sample[,col]/prior_std_list[col]
    }
  }
  ###########################################
  
  # For each number of bins (for each summary statistics S) - make ABC
  for (i in seq(1, n_bins_max)) {
    print(c("Begining to check sum.statistics with ", i, " n_bins......"))
    MRSSE <- 0.
    
    # For new data generated by each vector of nrow(sample) (n_obs) parameters
    # make ABC and find RSSE
    for (j in seq(1, nrow(sample))) {
      print(c("Using ", j, " of ", nrow(sample), " simulated parameter vectors"))
      par <- params[j,]
      new_data <- create_survey(data, mu_logv0=par[1], std_logv0=par[2],
                                mu_logtb=par[3], std_logtb=par[4], alpha_e=5,
                                beta_e=par[5])
    
      bins <- calc_bins(new_data, i)
      cur_model <- model2(new_data, bins)
      sum_stat_observed <- calc_fractions_in_bins(new_data, bins)
      ABC_rej_inner<-ABC_rejection(model=cur_model, prior=prior, nb_simul=nb_simul_MRSSE,
                                  summary_stat_target=sum_stat_observed, tol=tol_MRSSE,
                                  progress_bar=TRUE)
      abc_out_inner<-abc(sum_stat_observed, ABC_rej_inner$param,
                         ABC_rej_inner$stats, tol=0.5, method="neuralnet",
                         numnet=50, sizenet=2)
    
      # Make or not correction ##################
      if (use_correction==FALSE) {
       new_sample = ABC_rej_inner$param
      }
      else if (use_correction==TRUE) {
       new_sample = abc_out_inner$adj.values
      }
      ###########################################
  
      # TODO: Use scale(new_sample, scale=as.vector(prior_std_list))
      # Standardize #############################
      # if (is.null(standard)) {
      #   for (col in seq(1, ncol(new_sample))) {
      #     new_sample[,col]<-(new_sample[,col]-mean(new_sample[,col]))/sd(new_sample[,col])
      #   } 
      # }
      if (standard==TRUE) {
        # Standardize a-la Ewan et al. 
        # prior_std_list - list of std of priors for parameters
        # (0.15, 0.15, 6./sqrt(12.), 7.5/sqrt(12.), 29./sqrt(12.))
        for (col in seq(1, ncol(new_sample))) {
          new_sample[,col]<-new_sample[,col]/prior_std_list[col]
        }
      }
      ###########################################
      
      # Find RSSE ###############################
      # Here both new_sample and sample are standardized
      RSSE <- sqrt(sum((t(apply(new_sample, 1, '-', sample[j,])))**2) / nrow(new_sample))
      print(c("RSSE for ", j, " out of ", nrow(sample), " is ", RSSE))
      MRSSE <- MRSSE + RSSE
    }
    
    # Average RSSE for each of n_obs runs
    print(c("MRSSE for ", i, " out of ", n_bins_max, " is ", MRSSE / nrow(sample)))
    MRSSE_output[i] <- MRSSE / nrow(sample)
  }
  return(MRSSE_output)
}