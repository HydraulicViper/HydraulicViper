# Generate a Monte Carlo sample using Sobol' low-discrepancy quasi-random sequences
#
# makeMCSample
# Makes a Monte Carlo sample using Sobol' sequences
#
# Parameters:
#     n = number of samples to draw
#     vals = list describing distributions, each consisting of:
#         - a variable name 'var', 
#         - the distribution name, 'dist', e.g. 'unif' for Uniform (see R's distribution names)
#         - the distribution parameters, 'params' (names vary by distribution)
#
# Returns:
#   a data frame with an index column "n" and each distribution in a named column
#
# 
# install.packages("randtoolbox")
# install.packages("plyr")
library("randtoolbox")
library("plyr")

makeMCSample <- function(n, vals, p = 1) {
  # Packages to generate quasi-random sequences
  # and rearrange the data
  require(randtoolbox)
  require(plyr)
  
  # Generate a Sobol' sequence
  if(p ==2 ) {sob <- sobol(n, length(vals), seed = 1234, scrambling = 1)}else {sob <- sobol(n, length(vals), scrambling = 1)}
  
  # Fill a matrix with the values
  # inverted from uniform values to
  # distributions of choice
  samp <- matrix(rep(0,n*(length(vals)+1)), nrow=n)
  samp[,1] <- 1:n
  for (i in 1:length(vals)) {
    # i=1
    l <- vals[[i]]
    dist <- l$dist
    params <- l$params
    fname <- paste("q",dist,sep="")
    samp[,i+1] <- do.call(fname,c(list(p=sob[,i]),params))
  }
  
  # Convert matrix to data frame and add labels
  samp <- as.data.frame(samp)
  names(samp) <- c("n",laply(vals, function(l) l$var))
  return(samp)
}