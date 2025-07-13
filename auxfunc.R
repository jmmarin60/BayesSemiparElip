rescale.data <- function(x, robust = FALSE){
  #
  # Rescales a data matrix by subtracting (a robust estimate of) the mean and 
  # dividing each column by (a robust estimate of) the standard deviation.
  #
  if (robust == FALSE){
    mean.x <- apply(x, 2, mean)
    sd.x <- apply(x, 2, sd)
  }else{
    
  }
  x <- sweep(x, 2, mean.x)
  x <- sweep(x, 2, sd.x)
  return(list("x" = x, "mean.x" = mean.x, "sd.x" = sd.x))
}
#
set.grid <- function(x, npts = 100){
  #
  # Defines a grid for plotting.
  #
  if (is.vector(x)){
    x <- matrix(x, ncol = 1)
  }
  d <- ncol(x)
  grid <- matrix(rep(NA, (npts+1)*d), ncol = d)
  for (j in 1:d){
    minx <- min(x[ , j])
    maxx <- max(x[ , j])
    grid[, j] <- minx+c(0:npts)*(maxx-minx)/npts
  }
  return(grid)
}
#
make.symmetric <- function(M){
  #
  # Symmetrize a matrix.
  #
  return(0.5*(M+t(M)))
}
#
SALTsampler <- function(oldtheta,ssalt = 0.05){
  #
  # Internal function to generate a sample on the simplex.
  #
  k <- length(oldtheta)
  c <- sample(c(1:k),2)
  i <- c[1]
  l <- c[2]
  theta <- oldtheta
  theta[i] <- rnorm(1,log(oldtheta[i])-log(1-oldtheta[i]),ssalt)
  theta[i] <- ifelse(theta[i] < 0,exp(theta[i])/(1+exp(theta[i])),
                     1/(1+exp(-theta[i])))
  if (k == 2){
    theta[l] <- 1-theta[i]
  }else{
    theta[-i] <- (1-theta[i])*(oldtheta[-i]/(1-oldtheta[i]))
    theta[l] <- 1-sum(theta[-l])
  }
  logqratio <- log(theta[i])-log(oldtheta[i])+(k-1)*(log(1-theta[i])-
                                                       log(1-oldtheta[i]))
  return(list("theta" = theta,"logqratio" = logqratio))
}
#
calc.DIC <- function(logf, fhat, scale = 0){
  #
  # Calculates the DIC, expected deviance and pD given the estimated log-
  # likelihood (logf) and the estimated density of each datum (fhat).
  # If the data are rescaled, the scale (log of the determinant of the A matrix)
  # is included to give the DIC of the original data.
  #
  n <- length(fhat)
  DIC <- -4*(logf)+2*sum(log(fhat))+2*n*scale
  dev <- -2*(logf-n*scale)
  pD <-  DIC-dev
  return(list("DIC3" = DIC, "dev" = dev, "pD" = pD))
}
#
plot.fitted.density <- function(x, grid, f, kernel = FALSE, true = NULL, 
                      xlabels = "x", cols = c("red", "blue")){
  #
  # Overlays a histogram of x with a plot of the fitted density (and the kernel
  # density estimates and true density function).
  #
}
#
plot.bivar.density <- function(x, grid, f, true = NULL, labels = c("x", "y"),
                               p = c(0.01, 0.05, 0.5, 0.95, 0.99)){
  
}
