relliptical <- function(n,mu,Sigma,dist = "beta",parameters){
  #
  # Generates deviates from an elliptical distribution,
  # Y = mu + (A%*%S)*sqrt(W), where S is a spherical variable
  # and W follows the distribution defined in dist. A%*%t(A) = t(A)%*%A = Sigma.
  # Inputs:
  #   n = number of deviates to generate.
  #   mu = location parameter vector.
  #   Sigma = scale matrix.
  #   dist = distribution of W.  Available distributions:
  #     beta, exp, F, gamma, lnorm, weib.
  # parameters = parameters of the distribution defined in dist.
  #
  if (!require("MTS") == TRUE){
    install.packages("MTS")
  }
  A <- MTS::msqrt(Sigma)$mtxsqrt
  k <- length(mu)
  y <- matrix(rnorm(n*k),nrow=n)
  sumy2 <- apply(y^2,1,sum)
  y <- sweep(y,1,sqrt(sumy2),"/")
  y <- y%*%A
  if (dist == "beta"){
    g <- rbeta(n,parameters$shape1,parameters$shape2)
  }else if (dist == "exp"){
    g <- rexp(n,rate = parameters$rate)
  }else if (dist == "F"){
    g <- rf(n,parameters$df1,parameters$df2)
  }else if (dist == "gamma"){
    g <- rgamma(n,shape = parameters$shape,rate = parameters$rate)
  }else if (dist == "lnorm"){
    g <- rlnorm(n,parameters$mu,parameters$sigma)
  }else if (dist == "weib"){
    g <- rweibull(n,parameters$shape,parameters$scale)
  }
  y <- sweep(y,1,sqrt(g),"*")
  y <- sweep(y,2,mu,"+")
  return(y)
}



delliptical <- function(x, mu, A, distr = "Gamma", parameters, log = FALSE){
  #
  # Calculates the density function of elliptical data, x, where: 
  # x = mu + A S R, and S is a hyperspherical variable and R is a non-negative
  # variable.
  #
  # INPUTS
  # x = a vector or matrix of data.
  # mu = location vector.
  # A = scale matrix (or number).
  # distr = distribution of the variable R.
  # parameters = parameters of distr.
  # log = TRUE if the log density is required and FALSE otherwise.
  #
  # OUTPUT
  # f = the density or log density of x.
  #
  dim <- length(mu)
  if (is.vector(x)){x <- matrix(x, ncol = dim)}
  if (is.vector(A)){A <- matrix(A, ncol = 1)}
  if (distr == "Constant"){
    f <- dhyperspherical(sweep(x, 2, mu)%*%solve(A), dim = dim, 
                         r = parameters[1], log = TRUE)-log(det(A))
  }else{
    Sigma <- t(A)%*%A
    r <- sqrt(mahalanobis(x, mu, Sigma))
    f <- lgamma(dim/2)-log(2)-0.5*dim*log(pi)-log(det(A))+(1-dim)*log(r)
    if (distr == "Beta"){
      f <- f+dbeta(r, parameters[1], parameters[2], log = TRUE)
    }else if (distr == "Chisquare"){
      f <- f+dchisq(r, parameters[1], log = TRUE)
    }else if (distr == "Exponential"){
      f <- f+dexp(r, parameters[1], log = TRUE)
    }else if (distr == "F"){
      f <- f+df(r, parameters[1], parameters[2], log = TRUE)
    }else if (distr == "Gamma"){
      f <- f+dgamma(r, parameters[1], parameters[2], log = TRUE)
    }else if (distr == "Inverse.Gamma"){
      f <- f+dgamma(1/r, parameters[1], parameters[2], log = TRUE)-2*log(r)
    }else if (distr == "Lognormal"){
      f <- f+dlnorm(r, parameters[1], parameters[2], log = TRUE)
    }else if (distr == "Sqrt.Gamma"){
      f <- f+dgamma(r^2, parameters[1], parameters[2], log = TRUE)+log(2)+log(r)
    }else if (distr == "Uniform"){
      f <- f+dunif(r, parameters[1], parameters[2], log = TRUE)
    }else if (distr == "Weibull"){
      f <- f+dweibull(r, parameters[1], parameters[2], log = TRUE)
    }
  }
  if (log == FALSE){
    f <- exp(f)
  }
  return(f)
}


delliptical.marg <- function(orig.dim, x, mu, A, distr = "Gamma", parameters, 
                             m = 100, log = FALSE){
  #
  # Calculates the marginal density in an elliptical model.
  #
  dim <- length(mu)
  if (orig.dim < dim){
    stop("orig.dim must be >= length(mu)")
  }else if (orig.dim == dim){
    f <- delliptical(x, mu, A, distr, parameters, log = FALSE)
  }else{
    sqrt.b <- sqrt(rbeta(m, dim/2, (orig.dim-dim)/2))
    if (is.matrix(x)){
      n <- nrow(x)
    }else{
      n <- ifelse(dim == 1, length(x), 1)
    }
    f <- rep(0, n)
    for (i in 1:m){
      f <- f+delliptical(x, mu, A*sqrt.b[i], distr, parameters, log = FALSE)
    }
    f <- f/m
  }
  if (log == TRUE){
    f <- log(f)
  }
  return(f)
}

