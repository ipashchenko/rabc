library(gdata)
col.names <- c('source', 'bl', 's_thr', 'status')
setwd("/home/ilya/code/rabc/rabc")
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
    fractions[i] <- as.double(length(y$status))/length(statuses$status)
  }
  return(fractions)
}

#def flux_(b, v0, tb):
#  """
#    Flux of circular gaussian source with full flux ``v0`` and brightness
#    temperature ``tb`` at baseline ``b``.
#    :param b:
#        Baseline [baseline, ED]
#    :param v0:
#        Amplitude of component [Jy].
#    :param tb:
#        Brightness temperature of source [K].
#    :return:
#        Value of correlated flux.
#    """
#    b = b * 12742. * 10. ** 3
#    k = 1.38 * 10 ** (-23)
#    return v0 * np.exp(-math.pi * b ** 2. * v0 * 10 ** (-26) / (2. * k * tb))

flux<-function(x){
  x[1] = x[1] * 12742. * 10. ** 3
  k = 1.38 * 10 ** (-23)
  pi = 3.14159
  result = x[2] * exp(-pi * x[1] ** 2. * x[2] * 10 ** (-26) / (2. * k * x[3]))
  return(result)
}

# function must accept vector of model parameters and return summary statistics
toy_model<-function(x){
  
}