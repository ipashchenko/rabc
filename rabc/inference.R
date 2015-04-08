# Do ABC for geiven band, priors and model type
inference <- function(band, n_bins, path_to_data='/home/ilya/code/rabc/rabc/',
                      prior=list(c("normal", -0.43, 0.15),
                                 c("normal", 0.94, 0.15),
                                 c("unif", 7., 16.),
                                 c("unif", 0, 3.5),
                                 c("unif", 1, 30)),
                      nb_simul=50000, tol=0.05) {
  
  data_lookup = c("c"="data_to_R_Cband.txt", "l"="data_to_R_Lband.txt",
                  "k"="data_to_R_Kband.txt")
  source("/home/ilya/code/rabc/rabc/data_load.R")
  data = load_data(data_lookup[band], path=path_to_data)
  bins = calc_bins(data, n_bins)
  sum_stat_obs <- calc_fractions_in_bins(data, bins)
  
  source("/home/ilya/code/rabc/rabc/models.R")
  model <- model2(data, bins)
  
  library("EasyABC")
  library("abc")
  ABC_rej<-ABC_rejection(model=model, prior=prior, nb_simul=nb_simul,
                         summary_stat_target=sum_stat_obs, tol=tol,
                         progress_bar=TRUE)
  abc_out<-abc(sum_stat_obs, ABC_rej$param, ABC_rej$stats, tol=0.5,
               method="loclinear")
  return(c(ABC_rej, abc_out))
}