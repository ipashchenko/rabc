# Find kNN-distance for given sample
distances_k <- function(sample, k, standard=NULL, prior_std_list=NULL) {
  
  # Now msample[1,] - first element of sample
  msample <- matrix(sample, nrow=dim(sample)[1], ncol=dim(sample)[2])
  
  # Standardize
  if (is.null(standard)) {
    for (col in seq(1, ncol(msample))) {
      msample[,col]<-(msample[,col]-mean(msample[,col]))/sd(msample[,col])
    } 
  }
  else if (standard=='eaton') {
    # Standardize a-la Ewan et al. 
    # prior_std_list - list of std of priors for parameters
    # (0.15, 0.15, 6./sqrt(12.), 7.5/sqrt(12.), 29./sqrt(12.))
    for (col in seq(1, ncol(msample))) {
      msample[,col]<-msample[,col]/prior_std_list[col]
    }
  }
  
  distances <- vector(length = dim(sample)[1])
  for (i in seq(1, dim(sample)[1])) {
    element <- msample[i,]
    # Now melement[j,] - is i-th element for any j 
    melement <- t(matrix(rep(element, dim(msample)[1]), nrow=dim(sample)[2], ncol=dim(sample)[1]))
    diffs_sqrd <- apply(msample - melement, 2, function(x) x**2)
    sums_diffs_sqrd <- apply(diffs_sqrd, 1, sum)
    sqrt_sums_diffs_sqrd <- sqrt(sums_diffs_sqrd)
    # indx <- which(sqrt_sums_diffs_sqrd == sort(sqrt_sums_diffs_sqrd)[k])
    distances[i]<-sort(sqrt_sums_diffs_sqrd)[k] 
  }
  return(distances)
}

# Calculate entropy of sample using kNN-distances
entropy <- function(sample, n_par, k=4, prior_std_list=NULL) {
  # vector of k-th nearest neighbors for sample
  if (is.null(prior_std_list)) {
    distances <- distances_k(sample, k)
  }
  else {
    distances <- distances_k(sample, k, standard="eaton", prior_std_list=prior_std_list)
  }
  n <- dim(sample)[1]
  print(c("n", n))
  return(log(pi**(n_par/2)/gamma(n_par/2+1)) - digamma(k) + log(n) + (n_par/n)* sum(log(distances)))
}

