# Create distributions of ``Tb`` and ``v0``
band='c'
data_lookup = c("c"="data_to_R_Cband.txt", "l"="data_to_R_Lband.txt",
                "k"="data_to_R_Kband.txt")
source("/home/ilya/code/rabc/rabc/data_load.R")
data = load_data(data_lookup[band], path='/home/ilya/code/rabc/rabc/')
N=length(data[, 1])
components <- sample(1:2,prob=c(0.85,0.15),size=N,replace=TRUE)
mus <- c(12,13,5)
sds <- c(0.3,0.3)
tb_samples <- rnorm(n=N,mean=mus[components],sd=sds[components])
hist(tb_samples, breaks = 30)

create_tb <- function(prob, mus, sds, n) {
  components <- sample(1:2,prob=prob,size=n,replace=TRUE)
  tb_samples <- rnorm(n=n,mean=mus[components],sd=sds[components])
  return(tb_samples)
}

source("/home/ilya/code/rabc/rabc/simulate_survey.R")
source("/home/ilya/code/rabc/rabc/inference.R")
data_lookup = c("c"="data_to_R_Cband.txt", "l"="data_to_R_Lband.txt",
                "k"="data_to_R_Kband.txt")
source("/home/ilya/code/rabc/rabc/data_load.R")
band ='c'
data = load_data(data_lookup[band], path='/home/ilya/code/rabc/rabc/')
fake_data = create_survey(data, mu_logtb = c(12, 13.5), std_logtb = c(0.3, 0.3),
                          prob = c(0.85, 0.15))

sum_stat_obs <- 0.

source("/home/ilya/code/rabc/rabc/models.R")
model <- model22(fake_data)

library("EasyABC")
library("abc")
prior=list(c("normal", -0.43, 0.15),
         c("normal", 0.94, 0.15),
         c("unif", 7., 16.),
         c("unif", 0, 3.5),
         c("unif", 1, 30))
ABC_rej<-ABC_rejection(model=model, prior=prior, nb_simul=10000,
                       summary_stat_target=sum_stat_obs, tol=0.05,
                       progress_bar=TRUE)
# abc_out<-abc(sum_stat_obs, ABC_rej$param, ABC_rej$stats, tol=0.5,
#              method="neuralnet", numnet=50, sizenet=2)

