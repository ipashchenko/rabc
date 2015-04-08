library("EasyABC")
library("abc")
source("/home/ilya/code/rabc/rabc/enthropy.R")
source("/home/ilya/code/rabc/rabc/data_load.R")
source("/home/ilya/code/rabc/rabc/models.R")


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
                                 nb_simul=50000, tol=0.05) {
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
      ABC_rej<-ABC_rejection(model=cur_model, prior=prior, nb_simul=10000,
                             summary_stat_target=sum_stat_observed, tol=0.02,
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

#find_n_bins_MRSSE()