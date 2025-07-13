rni <- function(n, mu, A, distr = "Gamma", parameters){
  #
  # Generates random deviates from a normal/independent variable, x, such that: 
  # x = mu + A Z / sqrt(W), where Z is (multivariate) normal and W is a non-
  # negative variable.
  # 
  # INPUTS
  # n = number of data to generate.
  # mu = location vector.
  # A = scale matrix (or number).
  # distr = distribution of the variable W.
  # parameters = parameters of distr.
  #
  # OUTPUT
  # a list containing two variables: x and  w.
  # x = a matrix or vector of NI data.
  # w = a vector of non-negative data generated from distr.
  #
  if (distr == "Beta"){
    w <- rbeta(n, parameters[1], parameters[2])
  }else if (distr == "Chisquare"){
    w <- rchisq(n, parameters[1])
  }else if (distr == "Constant"){
    w <- rep(parameters[1], n)
  }else if (distr == "Exponential"){
    w <- rexp(n, parameters[1])
  }else if (distr == "F"){
    w <- rf(n, parameters[1], parameters[2])
  }else if (distr == "Gamma"){
    w <- rgamma(n, parameters[1], parameters[2])
  }else if (distr == "Inverse.Gamma"){
    w <- 1/rgamma(n, parameters[1], parameters[2])
  }else if (distr == "Lognormal"){
    w <- rlnorm(n, parameters[1], parameters[2])
  }else if (distr == "Sqrt.Gamma"){
    w <- sqrt(rgamma(n, parameters[1], parameters[2]))
  }else if (distr == "Uniform"){
    w <- runif(n, parameters[1], parameters[2])
  }else if (distr == "Weibull"){
    w <- rweibull(n, parameters[1], parameters[2])
  }
  dim <- length(mu)
  x <- matrix(rnorm(n*dim), ncol = dim)
  x <- sweep(sweep(x%*%A, 1, sqrt(w), "/"), 2, mu, "+")
  if (dim == 1 || n == 1){x <- as.vector(x)}
  return(list("x" = x, "w" = w))
}
#
# n <- 5
# mu <- 1
# A <- 2
# distr <- "Gamma"
# parameters <- c(1,2)
# rni(n, mu, A, distr, parameters)
# mu <- c(-1, 1)
# A <- matrix(c(3,1,1,2), nrow = 2)
# rni(n, mu, A, distr, parameters)
#
dni <- function(x, mu, A, distr = "Gamma", parameters, m = 100, log = FALSE){
  #
  # Calculates the density of x, where x is normal/independent:
  # x = mu + A Z / sqrt(W) where Z is (multivariate) normal and W is a non-
  # negative variable.
  # 
  # INPUTS
  # x = a vector or matrix of data.
  # mu = location vector.
  # A = scale matrix (or number).
  # distr = distribution of the variable W.
  # parameters = parameters of distr.
  # m = a Monte Carlo sample size for w.
  # log = TRUE if the log density is required and FALSE otherwise.
  #
  # OUTPUT
  # f = the density or log density of x.
  #
  dim <- length(mu)
  if (is.vector(x)){x <- matrix(x, ncol = dim)} 
  if (dim == 1){A <- matrix(A, ncol = 1)}
  Sigma <- t(A)%*%A
  if (distr == "Beta"){
    w <- rbeta(m, parameters[1], parameters[2])
  }else if (distr == "Chisquare"){
    w <- rchisq(m, parameters[1])
  }else if (distr == "Constant"){
    w <- rep(parameters[1], m)
  }else if (distr == "Exponential"){
    w <- rexp(m, parameters[1])
  }else if (distr == "F"){
    w <- rf(m, parameters[1], parameters[2])
  }else if (distr == "Gamma"){
    w <- rgamma(m, parameters[1], parameters[2])
  }else if (distr == "Inverse.Gamma"){
    w <- 1/rgamma(m, parameters[1], parameters[2])
  }else if (distr == "Lognormal"){
    w <- rlnorm(m, parameters[1], parameters[2])
  }else if (distr == "Sqrt.Gamma"){
    w <- sqrt(rgamma(m, parameters[1], parameters[2]))
  }else if (distr == "Uniform"){
    w <- runif(m, parameters[1], parameters[2])
  }else if (distr == "Weibull"){
    w <- rweibull(m, parameters[1], parameters[2])
  }
  f <- rep(0, dim(x)[1])
  for (i in 1:m){
    f <- f+LaplacesDemon::dmvn(x, mu, Sigma/w[i], log = FALSE)
  }
  f <- f/m
  if (log == TRUE){
    f <- log(f)
  }
  return(f)
}
#
# n <- 1000
# mu <- 1
# A <- 2
# distr <- "Gamma"
# parameters <- c(1,2)
# x <- rni(n, mu, A, distr, parameters)$x
# grid <- min(x)+c(0:1000)*(max(x)-min(x))/1000
# f <- dni(grid, mu, A, distr, parameters, log = FALSE)
# hist(x, breaks = 40, freq = FALSE)
# lines(grid, f)
# mu <- c(-1, 1)
# A <- matrix(c(3,1,1,2), nrow = 2)
# x <- rni(n, mu, A, distr, parameters)$x
# grid1 <- min(x[,1])+c(0:1000)*(max(x[,1])-min(x[,1]))/1000
# grid2 <- min(x[,2])+c(0:1000)*(max(x[,2])-min(x[,2]))/1000
# grid <- cbind(rep(grid1, 1001), rep(grid2, each = 1001))
# f <- dni(grid, mu, A, distr, parameters, log = FALSE)
# greencols <- c('greenyellow', 'green1', 'green2', 'green3', 'green4', 'seagreen')
# plot(hdrcde::hdr.2d(x = x[,1], y = x[,2], prob = c(1,5,50,95,99), 
#     den=list(x = grid1, y = grid2, z = matrix(f, nrow = 1001))), show.points = TRUE, shadecols = greencols)

