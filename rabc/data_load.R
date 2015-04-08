# Read & Write objects example
# On odin
# saveRDS(ABC_rej3, 'ABC_rej3.rds')
# on calculon
# ABC_rej3 <- readRDS("ABC_rej3.rds")

load_data <- function(fname, path="/home/ilya/code/rabc/rabc") {
  library(gdata)
  col.names <- c('source', 'bl', 's_thr', 'status')
  setwd(path)
  data = read.table(fname, col.names=col.names)
  return(data)
}

# Calculate fractions for different number of bins
calc_bins <- function(data, n_bins) {
  bins <- vector(length=n_bins+1)
  d = (max(data$bl) - min(data$bl)) / n_bins
  for (i in seq(1, n_bins+1)) {
    bins[i] <- min(data$bl)+(i-1)*d
  }
  return(bins)
}

# Calculate detection fractions in given baseline bins
calc_fractions_in_bins <- function(data, bins) {
  fractions <- vector(length=length(bins)-1)
  for (i in seq(1, length(bins)-1)) {
    statuses<-subset(data, bl>bins[i] & bl<bins[i+1], status)
    y<-subset(statuses, statuses == 'y')
    fractions[i] <- as.double(length(y$status))/length(statuses$status)
  }
  return(fractions)
}