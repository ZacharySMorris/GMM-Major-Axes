##function to identify unique pairs from list of items##
unique.pairs <- function(v){
  # v = a vector of group names to be paired
  
  n <- length(v)
  n_v <- c(1:n)
  
  temp.v <- lapply(n_v, FUN = function(x) (x+1):n)
  unique.v <- temp.v[-n]
}
#