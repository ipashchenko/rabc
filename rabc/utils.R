library(gdata)
col.names <- c('bsls', 's_thrs', 'status')
bsls_sthr_status = read.table('data.txt', col.names=col.names)

# Model that given parameters vector ``x`` returns vector of summary statistic
model <- function(x){
  
}

borders <- c(2., 5., 10., 15., 20.)
#test commit windows
det_fractions_in_bsl_ranges <- function(bsls_sthr_status, borders){
  for (i in seq(1, length(borders)-1)) {
    statuses = subset(dat, baseline_ed>borders[i] & baseline_ed<borders[i+1], status)
    fractions[i] <- sum(statuses)/nrow(statuses)
  }
  return(fractions)
}

