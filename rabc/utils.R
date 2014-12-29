library(gdata)
col.names <- c('source', 'bl', 's_thr', 'status')
data = read.table('data_to_R.txt', col.names=col.names)

# Model that given parameters vector ``x`` returns vector of summary statistic
model <- function(x){
  
}

borders <- c(2., 5., 10., 15., 20.)
#test commit windows
det_fractions_in_bsl_ranges <- function(data, borders){
  fractions <- vector(length=length(borders)-1)
  for (i in seq(1, length(borders)-1)) {
    statuses = subset(data, bl>borders[i] & bl<borders[i+1], status)
    y = subset(statuses, statuses == 'y')
    print(y)
    fractions[i] <- as.double(length(y))/length(statuses)
  }
  return(fractions)
}

toy_model<-function(x){
  
}