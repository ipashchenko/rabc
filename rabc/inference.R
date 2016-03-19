# Do ABC for geiven band, priors and model type
inference <- function(band, path_to_data='/home/ilya/code/rabc/rabc/',
                      prior=list(c("normal", -0.43, 0.15),
                                 c("normal", 0.94, 0.15),
                                 c("unif", 7., 16.),
                                 c("unif", 0, 3.5),
                                 c("unif", 1, 30)),
                      nb_simul=10000, tol=0.05) {

  data_lookup = c("c"="data_to_R_Cband.txt", "l"="data_to_R_Lband.txt",
                  "k"="data_to_R_Kband.txt")
  source("/home/ilya/code/rabc/rabc/data_load.R")
  data = load_data(data_lookup[band], path=path_to_data)
  sum_stat_obs <- 0.

  source("/home/ilya/code/rabc/rabc/models.R")
  model <- model22(data)

  library("EasyABC")
  library("abc")
  ABC_rej<-ABC_rejection(model=model, prior=prior, nb_simul=nb_simul,
                         summary_stat_target=sum_stat_obs, tol=tol,
                         progress_bar=TRUE)
  #abc_out<-abc(sum_stat_obs, ABC_rej$param, ABC_rej$stats, tol=0.5,
  #             method="neuralnet", numnet=50, sizenet=2)
  return(ABC_rej)
}

