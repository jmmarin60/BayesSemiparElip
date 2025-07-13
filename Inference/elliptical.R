rm(list = ls())

path1 = "~/Codigo/"
path2 = "~/Codigo/inference/"
path3 = "~/Codigo/distributions/"

fittedmodel = "elliptical"
truemodel = "ellip"
numrep = 100 # 100 samples
n = 200 # n=50, n=500

burnin <-10000
iters <- 100000
totits <- burnin+iters


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


#
# Bayesian inference for a (multivariate) elliptical distribution: 
# x[i] = mu + A S[i, ]*w[i],
# where S[ , ] is (multivariate) spherically uniform, mu is a location vector, 
# A is a scale matrix such that t(A)%*%A = Sigma, where Sigma is a symmetric, 
# positive definite matrix such that all of the diagonal elements sum to 1 and 
# w[i] is a non-negative, such that:
#
# Priors
#
# w[i] | aw, bw ~ Gamma(aw, bw),
# aw, bw | F ~ F
# F | alpha, F0 ~ DP(alpha F0),
# F0(aw, bw) = Exp(aw | cw), Exp(bw | dw),
# cw ~ Gamma(acw, bcw),
# dw ~ Gamma(adw, bdw),
# alpha ~ Gamma(aalpha, balpha),
# mu ~ N(m, V),
# sigma2 = diag(Sigma) ~ Dirichlet(asigma2),
# Write Sigma = diag(sqrt(sigma2))%*%R%*%diag(sqrt(sigma2)).
# Set W = D%*%R%*%D where D is a diagonal matrix.
# W ~ inverse Wishart(d+1, diag(d)*(d+1)), where d is the dension.
#
likelihood.elliptical <- function(x, mu, A, aw, bw, log = TRUE){
  #
  # Calculates the conditional likelihood function of x | w.
  #
  if (!exists("delliptical")){
    source(paste(path3, "distr.Elliptical.R", sep=""))
  }
  f <- sum(delliptical(x, mu, A, "Gamma", c(aw, bw), log = TRUE))
  if (log == FALSE){f <- exp(f)}
  return(f)
}
#
priors.elliptical <- function(d){
  #
  # Fixes priors.
  #
  priors <- list()
  priors$acw <- 1#0.01
  priors$bcw <- 1#0.01
  priors$adw <- 1#0.01
  priors$bdw <- 1#0.01
  priors$aalpha <- 2
  priors$balpha <- 2
  priors$m <- rep(0, d)
  priors$V <- diag(d)*10
  if (d > 1){
    priors$asigma2 <- 1
  }
  return(priors)
}
#
summary.stats.elliptical <- function(x){
  #
  # Calculates some summary statistics.
  #
  if (is.vector(x)){
    x <- matrix(x, ncol = 1)
  }
  return(list("n" = dim(x)[1], "sumx" = apply(x, 2, sum)))
}
#
calc.w.elliptical <- function(x, mu, Sigma){
  #
  # Calculates the value of w.
  #
  d <- length(mu)
  if (is.vector(x)){x <- matrix(x, ncol = d)}
  if (is.vector(Sigma)){Sigma <- matrix(Sigma, ncol = 1)}
  return(sqrt(mahalanobis(x, mu, Sigma)))
}
#
inits.elliptical <- function(x){
  #
  # Fixes initial values. 
  #
  inits <- list()
  if (is.vector(x)){x <- matrix(x, ncol = 1)}
  d <- ncol(x)
  inits$mu <- apply(x, 2, mean)
  #inits$mu <- rep(0,3)
  inits$Sigma <- cov(x)
  inits$Sigma <- inits$Sigma/sum(diag(inits$Sigma))
  inits$w <- calc.w.elliptical(x, inits$mu, inits$Sigma)
  b <- mean(inits$w)/var(inits$w)
  a <- b*mean(inits$w)
  inits$aw <- a*runif(n, 0.95, 1.05)
  inits$bw <- b*runif(n, 0.95, 1.05)
  inits$cw <- 0.1
  inits$dw <- 0.1
  inits$alpha <- 2
  if (d > 1){
    inits$sigma2 <- diag(inits$Sigma)
    inits$R <- cov2cor(inits$Sigma)
    inits$D <- diag(d)
    inits$W <- inits$D%*%inits$R%*%inits$D
  }
  return(inits)
}
#
generate.awbw.elliptical <- function(w, oldaw, oldbw, cw, dw, alpha){
  #
  # Generates from the conditional posterior distribution of aw, bw | w.
  #
  n <- length(w)
  aw <- oldaw
  bw <- oldbw
  for (i in 1:n){
    p <- dgamma(w[i], aw, bw)
    p[i] <- exp(log(alpha)+log(cw)+log(dw)-log(w[i])-log(w[i]+dw)-
                  2*log(cw-log(w[i])+log(w[i]+dw)))
    indic <- sample(c(1:n), 1, FALSE, p)
    if (indic == i){
      aw[i] <- rgamma(1, 2, cw-log(w[i])+log(w[i]+dw))
      bw[i] <- rgamma(1, aw[i]+1, w[i]+dw)
    }else{
      aw[i] <- aw[indic]
      bw[i] <- bw[indic]
    }
  }
  return(list("aw" = aw, "bw" = bw))
}
#
generate.cwdw.elliptical <- function(aw, bw, acw, bcw, adw, bdw){
  #
  # Generates from the conditional posterior distribution of cw | aw, bw | dw.
  #
  n.clust <- length(unique(aw))
  sum.aw <- sum(unique(aw))
  sum.bw <- sum(unique(bw))
  return(list("cw" = rgamma(1, acw+n.clust, bcw+sum.aw), 
              "dw" = rgamma(1, adw+n.clust, bdw+sum.bw), "n.clust" = n.clust))
}
#
generate.alpha.elliptical <- function(oldalpha, aalpha, balpha, n, n.clust){
  #
  # Generates from the conditional posterior of alpha.
  #
  theta <- rbeta(1, oldalpha+1, n)
  c <- aalpha+n.clust-1
  d <- balpha-log(theta)
  p <- c(c, n*d)
  alpha <- sample(c(rgamma(1, c+1, d), rgamma(1, c, d)), size = 1, 
                  replace = FALSE, prob = p)
  return(alpha)
}
#
# generate.mu.elliptical <- function(x, oldmu, Sigma, aw, bw, oldloglik, m,
#                                    invV, par = 0.0001){
#   #
#   # Generates from the conditional posterior distribution of mu | x, Sigma, w.
#   #
#   acc <- 1
#   mu <- rnorm(length(oldmu), oldmu, par)
#   A <- amen::mhalf(Sigma)
#   loglik <- likelihood.elliptical(x, mu, A, aw, bw, TRUE)
#   paccept <- loglik-oldloglik+LaplacesDemon::dmvn(mu, m, V, TRUE)-
#     LaplacesDemon::dmvn(oldmu, m, V, TRUE)
#   if (log(runif(1)) > paccept){
#     acc <- 0
#     mu <- oldmu
#     loglik <- oldloglik
#   }
#   return(list("mu" = mu, "loglik" = loglik, "acc" = acc))
# }

generate.mu.elliptical <- function(x, mu, Sigma, aw, bw, loglik, par.mu = 0.01){
  mu.sd = rep(par.mu,length(mu))
  V = diag(length(mu))
  m = rep(0,length(mu))
  if (!require("mvtnorm")){
    install.packages("mvtnorm")
  }
  #
  # Generates from the conditional posterior distribution of mu.
  #
  k <- length(mu)
  old.loglik <- loglik
  old.mu <- mu
  acc <- rep(1, k)
  for (i in 1:k){
    mu[i] <- rnorm(1, old.mu[i], mu.sd[i])
    log.paccept <- likelihood.elliptical(x, mu, amen::mhalf(Sigma), aw, bw, TRUE) - old.loglik + ifelse(k == 1, dnorm(mu, m, sqrt(V), TRUE)-
dnorm(old.mu, m, sqrt(V), TRUE),mvtnorm::dmvnorm(mu, m, V, TRUE)-mvtnorm::dmvnorm(old.mu, m, V, TRUE))
    if (log(runif(1)) > log.paccept){
      loglik <- old.loglik
      mu <- old.mu
      acc[i] <- 0
    }else{
      old.loglik <- loglik
      old.mu <- mu
    }
  }
  return(list("loglik" = loglik, "mu" = mu, "acc" = acc))
}

#
generate.sigma2.elliptical <- function(x, mu, oldSigma, R, aw, bw, oldloglik,  
                                       asigma2, par = 0.025){
  #
  # Generates from the conditional posterior distribution of sigma2 | x, mu, R, 
  # w.
  #
  acc <- 1
  oldsigma2 <- diag(oldSigma)
  result <- SALTsampler(oldsigma2, par)
  sigma2 <-result$theta
  Sigma <- diag(sqrt(sigma2))
  Sigma <- Sigma%*%R%*%Sigma
  A <- amen::mhalf(Sigma)
  loglik <- likelihood.elliptical(x, mu, A, aw, bw, TRUE)
  paccept <- loglik-oldloglik+result$logqratio+
    (asigma2-1)*sum(log(sigma2)-log(oldsigma2))
  if (log(runif(1)) > paccept){
    acc <- 0
    Sigma <- oldSigma
    loglik <- oldloglik
  }
  return(list("Sigma" = Sigma, "loglik" = loglik, "acc" = acc))
}
#
generate.R.elliptical <- function(x, mu, oldSigma, oldR, oldD, oldW, aw, bw,  
                                  oldloglik, par = 500){
  #
  # Generates from the conditional posterior distribution of R | x, mu, sigma2, 
  # w.
  #
  acc <- 1
  W <- LaplacesDemon::rwishart(par,oldW/par)
  D <- diag(sqrt(diag(W)))
  R <- cov2cor(W)
  sigma <- diag(sqrt(diag(oldSigma)))
  Sigma <- sigma%*%R%*%sigma
  A <- amen::mhalf(Sigma)
  loglik <- likelihood.elliptical(x, mu, A, aw, bw, TRUE)
  paccept <- loglik-oldloglik+
    LaplacesDemon::dinvwishart(W, d+1, diag(d)*(d+1), log = TRUE)-
    LaplacesDemon::dinvwishart(oldW, d+1, diag(d)*(d+1), log = TRUE)+
    LaplacesDemon::dwishart(oldW, par, W/par,log = TRUE)-
    LaplacesDemon::dwishart(W, par,oldW/par,log = TRUE) 
  if (log(runif(1)) > paccept){
    acc <- 0
    Sigma <- oldSigma
    R <- oldR
    D <- oldD
    W <- oldW
    loglik <- oldloglik
  }
  return(list("Sigma" = Sigma, "R" = R, "D" = D, "W" = W, "loglik" = loglik, 
              "acc" = acc))
}
#
integrated.delliptical.marg <- function(x, mu, A, caw, daw, cbw, dbw){
  #
  # Calculates the marginal elliptical distribution with the parameters aw, bw
  # integrated out.
  #
}
#
# SAMPLE CODE.
#
# Load required functions.
#
source(paste(path1, "auxfunc.R", sep=""))
source(paste(path3, "distr.Elliptical.R", sep=""))
# ------------------------------------------------------------------------------
#
# Fix settings for scaling, graphs, DIC and MCMC.
#
simulate.data <- TRUE
plots <- FALSE
if (plots == TRUE){npts <- 100}
dic <- TRUE

# setwd(paste(path1, "experiments", sep=""))
# dir.create(paste(fittedmodel,"_",truemodel,"_n", n, "nrep", numrep,sep=""))
# setwd(paste0(getwd(),paste("/",paste(fittedmodel,"_",truemodel,"_n", n, "nrep", numrep,sep=""),sep="")))
# getwd()
seed_vector = 1:numrep
DIC3 = vector(length = numrep+1)
DEV = vector(length = numrep+1)
pD = vector(length = numrep+1)
tiempomatriz = matrix(nrow = numrep+1, ncol = 3)

## Load data
for (jindex in 1:numrep){
  #print(jindex)
  set.seed(seed_vector[jindex])
  true.mu <- c(2,10,5)
  true.Sigma <- matrix(c(5,2,2,2,5,2,2,2,5), ncol = 3, byrow = TRUE)
  if (is.vector(true.Sigma)){true.Sigma <- matrix(true.Sigma, ncol = 1)}
  true.w <- sum(diag(true.Sigma))
  true.scaled.Sigma <- true.Sigma/true.w
  x <- relliptical(n=n,mu=true.mu,Sigma = true.Sigma, dist = "weib",parameters = list(shape = 5,scale = 2))
  
  
  # Load or generate data.
  #
  # if (simulate.data == TRUE){
  #   n <- 1000
  #   true.mu <- c(-3,5,10)
  #   true.Sigma <-  matrix(c(3,1,1,1,3,1,1,1,3), ncol = 3)
  #   if (is.vector(true.Sigma)){true.Sigma <- matrix(true.Sigma, ncol = 1)}
  #   if (!require("amen")){install.packages("amen")}
  #   true.A <- amen::mhalf(true.Sigma)
  #   true.distr <- "Gamma"
  #   true.pars <- c(1,2)
  #   source("C:/Research/2023-2024/Elliptically contoured distributions/R codes/Distributions/distr.Hyperspherical.R")
  #   source("C:/Research/2023-2024/Elliptically contoured distributions/R codes/Distributions/distr.Elliptical.R")
  #   if (!require("amen")){install.packages("amen")}
  #   result <- relliptical(n, true.mu, true.A, true.distr, true.pars)
  #   x <- result$x
  # }else{
  #   path <- "C:/Research/2023-2024/Elliptically contoured distributions/Data/data.csv"
  #   x <- read.csv(path, sep = ",")
  # }
  #
  
  
  
  
  # If x is a vector, make it a matrix and calculate the number of variables.
  #
  if (is.vector(x)){x <- matrix(x, ncol = 1)}
  n <- nrow(x)
  d <- ncol(x)
  #
  # Settings for MCMC.
  #
  par.awbw <- 0.025
  awbw.acc <- 0 
  par.mu <- 0.001
  mu.acc <- 0
  if (d > 1){
    par.sigma2 <- 0.001
    sigma2.acc <- 0
    par.R <- 500
    R.acc <- 0
  }
  simuls <- 10
  #
  # ------------------------------------------------------------------------------
  #
  # Set priors.
  #
  priors <- priors.elliptical(d)
  acw <- priors$acw
  bcw <- priors$bcw
  adw <- priors$adw
  bdw <- priors$bdw
  aalpha <- priors$aalpha
  balpha <- priors$balpha
  m <- priors$m
  V <- priors$V
  invV <- solve(V)
  if (d > 1){
    asigma2 <- priors$asigma2
  }
  rm(priors)
  #
  # Set initial values.
  #
  inits <- inits.elliptical(x)
  w <- inits$w
  aw <- inits$aw
  bw <- inits$bw
  n.clust <- length(unique(aw))
  clust.size <- rep(NA, n)
  cw <- inits$cw
  dw <- inits$dw
  alpha <- inits$alpha
  mu <- inits$mu
  Sigma <- inits$Sigma
  if (is.vector(Sigma)){Sigma <- matrix(Sigma, ncol = 1)}
  invSigma <- solve(Sigma)
  if (d > 1){
    sigma2 <- inits$sigma2
    R <- inits$R
    D <- inits$D
    W <- inits$W
    W <- make.symmetric(W)
    A <- amen::mhalf(Sigma)
  }else{
    A <- sqrt(Sigma[1, 1])
  }
  loglik <- likelihood.elliptical(x, mu, A, aw, bw, log = TRUE)
  #
  # Initialize posterior statistics.
  #
  ew <- 0
  emu <- rep(0, d)
  eSigma <- matrix(rep(0, d^2), ncol = d)
  #
  
  mupost = rep(NA, d)
  wpost = c()
  Sigma11post = c()
  Sigma22post = c()
  Sigma33post = c()
  Sigma12post = c()
  Sigma13post = c()
  Sigma23post = c()
  alphapost=c()
  if (plots == TRUE){
    #
    # Set grid for plots, calculate the true density (if the data are simulated) 
    # and initialize the predictive densities.
    #
    npts.plus.1 <- npts+1
    grid <- set.grid(x, npts)
    xlabels <- c(expression(x[1]), expression(x[2]), expression(x[3]), 
                 expression(x[4]), expression(x[5]), expression(x[6]))
    line.col <- "green"
    marg.f <- matrix(rep(0, npts.plus.1*d), ncol = d)
    if (d > 1){
      if (!require("hdrcde") == TRUE){
        install.packages("hdr.cde")
      }
      contour.cols <- c('greenyellow','green1','green2','green3','green4',
                        'seagreen')
      d.minus.1 <- d-1
      bivar.f <- array(rep(0, (d*(npts.plus.1))^2),  dim = c(npts.plus.1, 
                                                             npts.plus.1, d, d))
    }
    if (simulate.data == TRUE){
      #
      # Calculate the true (marginal and bivariate) densities.
      #
      true.line.col <- "red"
      true.marg.f <- matrix(nrow = npts.plus.1, ncol = d)
      for (j in 1:d){
        true.marg.f[ , j] <- delliptical.marg(d, grid[ , j], true.mu[j],  
                                              true.A[j, j], true.distr, true.pars, 100)
      }
      if (d > 1){
        true.bivar.f <- array(dim = c(npts.plus.1, npts.plus.1, d, d))
        for (i in 1:d.minus.1){
          for (j in (i+1):d){
            true.bivar.f[, , i, j] <- 
              delliptical.marg(d, cbind(rep(grid[,i], npts.plus.1),
                                        rep(grid[,j], each = npts.plus.1)), c(true.mu[i], true.mu[j]), 
                               amen::mhalf(matrix(c(true.Sigma[i, i], true.Sigma[i, j], 
                                                    true.Sigma[i, j], true.Sigma[j, j]), nrow = 2)), 
                               true.distr, true.pars, 100)
          }
        }
      }
    }
  }
  #
  if (dic == TRUE){
    dev <- 0
    fhat <- rep(0, n)
    # if (simulate.data == TRUE){
    #   true.dev <- -2*likelihood.elliptical(x, true.mu, true.A, true.pars[1], 
    #                                        true.pars[2], TRUE)
    # }
  }
  #
  # Start sampler.
  #
  for (i in 1:totits){
    i
    w <- calc.w.elliptical(x, mu, Sigma)
    result <- generate.awbw.elliptical(w, aw, bw, cw, dw, alpha)
    aw <- result$aw
    bw <- result$bw
    loglik = likelihood.elliptical(x, mu, amen::mhalf(Sigma), aw, bw, TRUE)
    result <- generate.cwdw.elliptical(aw, bw, acw, bcw, adw, bdw)
    cw <- result$cw
    dw <- result$dw
    n.clust <- result$n.clust
    alpha <- generate.alpha.elliptical(alpha, aalpha, balpha, n, n.clust)
    result <- generate.mu.elliptical(x, mu, Sigma, aw, bw, loglik, par.mu)
    mu <- result$mu
    loglik <- result$loglik
    mu.acc <- mu.acc+result$acc
    if (d > 1){
      result <- generate.sigma2.elliptical(x, mu, Sigma, R, aw, bw, loglik,  
                                           asigma2, par.sigma2)
      Sigma <- result$Sigma
      loglik <- result$loglik
      sigma2.acc <- sigma2.acc+result$acc
      result <- generate.R.elliptical(x, mu, Sigma, R, D, W, aw, bw, loglik, 
                                      par.R)
      Sigma <- result$Sigma
      invSigma <- solve(Sigma)
      R <- result$R
      D <- result$D
      W <- result$W
      loglik <- result$loglik
      R.acc <- R.acc+result$acc
    }
    if (i > burnin){
      ew <- ew+rgamma(1, aw, bw)
      emu <- emu+mu
      eSigma <- eSigma+Sigma
      unique.aw <- unique(aw)
      unique.bw <- unique(bw)
      alpha.plus.n <- alpha+n
      for (l in 1:n.clust){
        clust.size[l] <- sum(aw == unique.aw[l])
      }
      pred.aw <- rexp(simuls, cw)
      pred.bw <- rexp(simuls, dw)
      if (plots == TRUE){
        for (j in 1:d){
          result <- 0
          for (l in 1:n.clust){
            result <- result+clust.size[l]*delliptical.marg(d, grid[ , j], mu[j], 
                                                            sqrt(Sigma[j, j]), "Gamma", c(aw[l],bw[l]), m = 10) 
          }
          marg.f[ , j] <- marg.f[ , j]+result/alpha.plus.n
          result <- 0
          for (m in 1:simuls){
            result <- result+delliptical.marg(d, grid[ , j], mu[j], 
                                              sqrt(Sigma[j, j]), "Gamma", c(pred.aw[m],pred.bw[m]), m = 10)
          }
          marg.f[ , j] <- marg.f[ , j]+alpha*result/(simuls*alpha.plus.n)
        }
        if (d > 1){
          for (j in 1:d.minus.1){
            for (k in (j+1):d){
              S <- amen::mhalf(matrix(c(Sigma[j, j], Sigma[j, k], Sigma[j, k], 
                                        Sigma[k, k]), nrow = 2))
              result <- 0
              for (l in 1:n.clust){
                result <- result+clust.size[l]*delliptical.marg(d, 
                                                                cbind(rep(grid[ , j], npts.plus.1), rep(grid[ , k], 
                                                                                                        each = npts.plus.1)), c(mu[j], mu[k]), S, "Gamma", 
                                                                c(aw[l], bw[l]), m = 10)  
              }
              bivar.f[ , , j, k] <- bivar.f[ , , j, k]+result/alpha.plus.n
              result <- 0
              for (m in 1:simuls){
                result <- result+delliptical.marg(d, cbind(rep(grid[ , j], 
                                                               npts.plus.1), rep(grid[ , k], each = npts.plus.1)), 
                                                  c(mu[j], mu[k]), S, "Gamma", c(pred.aw[m], pred.bw[m]), m = 10)  
              }
              bivar.f[ , , j, k] <- bivar.f[ , , j, k]+alpha*result/
                (simuls*alpha.plus.n)
            }
          }
        }
      }
      if (dic == TRUE){
        if (d == 1){
          result <- 0
          result1 <- 0
          for (l in 1:n.clust){
            result <- result+clust.size[l]*delliptical(x, mu, sqrt(Sigma[1, 1]), 
                                                       "Gamma", c(aw[l], bw[l]))
          }
          gx = mahalanobis(x, mu, Sigma[1,1])
          result1 = log(cw) + log(dw)-log(sqrt(gx))- log(sqrt(gx) + dw) - 2*log(cw-log(sqrt(gx))+log(sqrt(gx)+dw))
          result1 = exp(result1)
          result <- (result+alpha*result1)/alpha.plus.n
        }else{
          Sigma <- 0.5*(Sigma+t(Sigma))
          A <- amen::mhalf(Sigma)
          result <- 0
          result1 <- 0
          for (l in 1:n.clust){
            result <- result+clust.size[l]*delliptical(x, mu, A, "Gamma", 
                                                       c(aw[l], bw[l]))
          }
          gx = mahalanobis(x, mu, A)
          result1 = log(cw) + log(dw)-log(sqrt(gx))- log(sqrt(gx) + dw) - 2*log(cw-log(sqrt(gx))+log(sqrt(gx)+dw))
          result1 = exp(result1)
          result <- (result+alpha*result1)/alpha.plus.n
        }
        dev <- dev+sum(log(result))
        fhat <- fhat+result
      }
      mupost = cbind(mupost,mu)
      alphapost = c(alphapost, alpha)
      Sigma11post = c(Sigma11post, Sigma[1,1])
      Sigma22post = c(Sigma22post, Sigma[2,2])
      Sigma33post = c(Sigma33post, Sigma[3,3])
      Sigma12post = c(Sigma12post, Sigma[1,2])
      Sigma13post = c(Sigma13post, Sigma[1,3])
      Sigma23post = c(Sigma23post, Sigma[2,3])
    }
  }
  
  #
  # Normalize output.
  # 
  ew <- ew/iters
  emu <- emu/iters
  eSigma <- eSigma/iters
  #awbw.acc <- awbw.acc/totits
  if (d > 1){
    sigma2.acc <- sigma2.acc/totits
    R.acc <- R.acc/totits
  }
  
  # mupost = mupost[,-1]
  # save(mupost, file = "mupost.RData")
  # Sigmapost = cbind(Sigma11post, Sigma22post, Sigma33post, Sigma12post, Sigma13post, Sigma23post)
  # save(Sigmapost, file = "Sigmapost.RData")
  # save(alphapost, file="alphapost.RData")
  # save(wpost, file="wpost.RData")
  # save(ew, file="ew.RData")
  # save(emu, file="emu.RData")
  # save(eSigma, file = "eSigma.RData")  
  
  
  #
  if (plots == TRUE){
    #
    # Draw plots.
    #
    marg.f <- marg.f/iters
    for (j in 1:d){
      hist(x[ , j], freq = FALSE, breaks = 20, xlab = xlabels[j], ylab = "f", 
           main = "")
      lines(grid[ , j], marg.f[ , j], lwd = 2, col = line.col)
      if (simulate.data == TRUE){
        lines(grid[ , j], true.marg.f[ , j], lwd = 2, col = true.line.col)
      }
    }
    if (d > 1){
      bivar.f <- bivar.f/iters
      for (j in 1:d.minus.1){
        for (k in (j+1):d){
          plot(hdrcde::hdr.2d(x = x[ , j], y = x[ , k], prob=c(1,5,50,95,99),
                              den=list(x = grid[ , j], y = grid[ , k], z = bivar.f[ , , j, k])),
               show.points = TRUE, xlab = xlabels[j], ylab = xlabels[k], 
               shadecols = contour.cols)
        }
      }
    }
  }
  #
  if (dic == TRUE){
    #
    # Calculate DIC3.
    #
    DEV[jindex] <- -2*dev/iters
    write.table(DEV, file = "dev.txt", sep =" ")
    fhat <- fhat/iters
    DIC3[jindex] <- 2*DEV[jindex]+2*sum(log(fhat))
    write.table(DIC3, file = "DIC3.txt", sep =" ")
    pD[jindex] <- DIC3[jindex]-DEV[jindex]
    write.table(pD, file = "pD.txt", sep =" ")
  }
}


mean_dic3 = mean(DIC3[1:jindex])
DIC3[jindex+1] = mean_dic3
#write.table(DIC3, file = "DIC3.txt", sep =" ")
mean_dev = mean(DEV[1:jindex])
DEV[jindex+1] = mean_dev
#write.table(DEV, file = "dev.txt", sep =" ")
mean_pD = mean(pD[1:jindex])
pD[jindex+1] = mean_pD
#write.table(pD, file = "pD.txt", sep =" ")

