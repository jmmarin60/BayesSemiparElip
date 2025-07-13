rhyperspherical <- function(n, dim = 1, r = 1){
  #
  # Generates data from a hyperspherical distribution.
  # 
  # INPUTS
  # n = number of data to generate.
  # dim = dimension of the hypersphere.
  # r = radius of the hypersphere.
  #
  # OUTPUT
  # x = a vector or matrix of hyperspherical data.
  #
  x <- matrix(rnorm(n*dim), nrow = n)
  x <- r*sweep(x, 1, apply(x, 1, norm, type = "2"), "/")
  if (dim == 1){x <- as.vector(x)}
  return(x)
}
# rhyperspherical(5, 1)
# (x <- rhyperspherical(5, 2, 2))
# The radius is 2 so the row sums of squares should all be 4.
# apply(x^2, 1, sum) 
#
dhyperspherical <- function(x, dim = ifelse(is.vector(x), 1, ncol(x)), r = 1, 
                            log = FALSE){
  #
  # Calculates the hyperspherical density function.
  #
  # INPUTS
  # x = a vector or matrix of data.
  # dim = dimension of the hypersphere.
  # r = radius of the hypersphere.
  #
  # OUTPUT
  # f = the density function of x.
  #
  if (is.vector(x)){x <- matrix(x, ncol = dim)}
  f <- rep(0, nrow(x))
  constraint <- abs(apply(x^2, 1, sum)-r^2)
  f[constraint < .Machine$double.eps^0.25] <- 0.5*gamma(dim/2)/((pi^(dim/2))*
                                                             (r^(dim-1)))
  if (log == TRUE){f <- log(f)}
  return(f)
}


# x <- rhyperspherical(5, 1)
# dhyperspherical(x, 1)
# x <- rhyperspherical(5, 2, 2)
# dhyperspherical(x, 2, 2)
# x[1, 1] <- 0
# dhyperspherical(x, 2, 2)
