
# Dependencies
require(CVXR) # This require using the MOSEK solver, follow instructions for downloading an academic license through https://cvxr.rbind.io/cvxr_examples/cvxr_using-other-solvers/
#require(ggplot2) # 


# Must be downloaded from https://github.com/richardkwo/multChernoff
require(multChernoff)



##   Basic math functions --------------------------------------------------

logit <- function(x){
  if ({
    any(x < 0)
  } || {
    any(x > 1)
  }) 
    stop("x must be in [0,1]")
  out = log(x/(1-x))
  return(out)
}

logistic <- function(x){
  return(1/(1 + exp(-x)))
}



##   Math functions  --------------------------------------------------
# Function: computes the kl divergence between any two distribution vectors p and q
# Input: p,q, both are probability distributions, can be numeric, matrix or array   
# Output: kl divergence between p and q 

kl_divergence <- function(p,q){
  if(!(all(dim(p) == dim(q)))){
    stop("p,q must be of the same dimensions" )
  }
  out <- p*log(p/q)
  out[is.nan(out)] <- 0
  return(sum(out))
}




# Function: computes the total variation distance between the empirical distributions of the samples 
# Input: x,y, two samples of a discrete distribution which need not be the same length,  numeric
# Output: total variation distance between the empirical distributions of x and y 
tv_norm <- function(x,y){
  n <- min(min(x), min(y))
  m <- max(max(x), max(y))
  
  Nx <- length(x)
  Ny <- length(y)
  
  series_term <- 0
  for(k in n:m){
    series_term <- series_term + abs((sum(x == k)/Nx) - (sum(y == k)/Ny))
  }
  out <- (1/2)*series_term
  return(out)
}


# Function: conditional distribution generated according to the measurement kernel model
# Input: N;  must be an integer 
#        ker, number of repeated observations per individual; integer 
#        h, bandwidth; > 0
# Output: Conditional Distribution Function for the Measurement Kernel Model
conditional_mkm <- function(N, ker, h){
  #ensuring these don't change when changing definitions of N,ker and h
  force(N)
  force(ker)
  force(h)
  i <- 0:N
  out.function <- function(gam){
    out <- ker((i - N*gam)/h)/sum(ker((i - N*gam)/h))
    return(out)
  }
  return(out.function)
}


# Function: compute empirical distribution function
# Input: N, (support size {0,...,N}); integer 
#        x, sample of observations; vector of integers
# Output: Empirical distribution function for x
compute_edf <- function(x,N){
  p.hat <- c()
  for(y in 0:(N)){
    prop <- mean(x == y)
    p.hat <- c(p.hat, prop)
  }
  return(p.hat)
}






# Function: computes a linear approximation to the cdf of a random variable, based a set of quantiles 
#           and corresponding inverse quantiles 
# Input: x: argument of the cdf; numeric 
#        quantiles:  vector of quantiles; numeric 
#        tau: vector of inverse quantiles; numeric 
# Output: Empirical distribution function for x
cdf_L <- function(x,quantiles, tau){
  
  i.set <- sapply(x, function(z){
    out <-  min(which(z <= quantiles))
    return(out)
  })
  
  out <- sapply(1:length(i.set), function(z){
    i = i.set[z]
    out <-  tau[i-1] + (tau[i] - tau[i-1])*((x[z] - quantiles[i-1])/(quantiles[i] - quantiles[i-1]))
    if(x[z] == 0 ){
      out <- 0
    } else if (x[z] == 1){
      out <- 1
    } else if (x[z]  < 0 | x[z]  > 1){
      out <- NA
    }
    return(out)
  })
  out <- unlist(out)
  length(out)
  return(out)
}

# Function: computes a linear approximation to the inverse cdf (or quantile function) 
#           of a random variable, based a set of quantiles and corresponding inverse quantiles 
# Input: x: argument of the cdf; numeric 
#        quantiles:  vector of quantiles; numeric 
#        tau: vector of inverse quantiles; numeric 
# Output: quantile function function for t
qf_L <- function(t,quantiles, tau){
  
  i.set <- sapply(t, function(z){
    out <-  min(which(z <= tau))
    return(out)
  })
  
  out <- sapply(1:length(i.set), function(z){
    i = i.set[z]
    out <-  quantiles[i-1] + (quantiles[i] - quantiles[i-1])*((t[z] - tau[i-1])/(tau[i] - tau[i-1]))
    if(t[z] == 0 ){
      out <- 0
    } else if (t[z] == 1){
      out <- 1
    } else if (t[z]  < 0 | t[z]  > 1){
      out <- NA
    }
    return(out)
  })
  
  
  
  
  return(out)
}



# Function: computed the implied distribution on the marginal from a latent distribution and a conditional
# Input: tau: vector of inverse quantiles; numeric 
#        latent.trait.quantiles:  vector of quantiles; numeric 
#        cond: function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        numeric.points: tuning parameter for the precision of the numerical approximation; integer 
# Output: p_ma: a probabilty vector on the observed scores 
compute_p_ma <- function(tau, latent.trait.quantiles,
                         cond, numeric.points = 100){
  
  R_bins <- length(latent.trait.quantiles) - 1
  
  A.matrix <- matrix(data = NA, nrow = N + 1, ncol = R_bins)
  for(i in 1:R_bins){
    design.points <- seq(latent.trait.quantiles[i], 
                         latent.trait.quantiles[i+1], 
                         length.out = numeric.points)
    
    A.row <- rep(0, N+1)
    for(j in 1:numeric.points){
      y = 0:N
      weights <- cond(design.points[j])
      A.row <- A.row + weights/(sum(weights))
    }
    A.row <- A.row/sum(A.row)
    A.matrix[,i] = A.row
  }
  weights <- c()
  for(i in 1:R_bins){
    weights[i] <- tau[i + 1] - tau[i]
  }
  p.ma <- A.matrix %*% weights
  return(p.ma)
}



# Function: computes the inverse quantiles from a set of uniformly binned latent weights 
# Input: latent.mixture: vector of weights in a uniform bin 
# Output: tau: a vector of the inverse quantiles corresponding to the edges of the latent mixture bins 

inv_quantiles_from_weights <- function(latent.mixture){
  tau <- rep(NA, length(latent.mixture) + 1)
  tau[1] <- 0
  run.sum <- 0
  for(i in 1:length(latent.mixture)){
    run.sum <- run.sum + latent.mixture[i]
    tau[i + 1] <- run.sum
  }
  # numerical error correction
  tau[tau > 1] = 1
  return(tau)
}


##   Kernel Functions ----------------------------------------

# all of these are simple functions which can be used to define a kernel function 

gaussian_kernel <- function(x){
  return(exp(-x^2/2))
}

exponential_kernel <- function(x){
  return(exp(-abs(x)))
}

logistic_kernel <- function(x){
  return(1/(exp(x) + 2 + exp(-x)))
}


## Kernel support  |x| <= 1
triangle_kernel <- function(x){
  out <- (1 - abs(x))*(abs(x) <= 1)
  return(out)
}


epanechnikov_kernel <- function(x){
  out <- (3/4)*(1 - x^2)*(abs(x) <= 1)
  return(out)
}

### Not Continuous which may cause problems 
uniform_kernel <- function(x){
  out <- 1*(abs(x) <= 1)
  return(out)
}

quartic_kernel <- function(x){
  out <- (15/16)*(1 - x^2)^2*(abs(x) <= 1)
  return(out)
}

triweight_kernel <- function(x){
  out <- (35/32)*(1 - x^2)^3*(abs(x) <= 1)
  return(out)
}

tricube_kernel <- function(x){
  out <- (70/81)*(1 - abs(x)^3)^3*(abs(x) <= 1)
  return(out)
}

cosine_kernel <- function(x){
  out <- (pi/4)*(cos(pi*x/2))^3*(abs(x) <= 1)
  return(out)
}




##   Latent model fitting functions  ----------------------------------------


# Function: sample from the latent distribution given an observed score. 
# Input: y.obs, the observed score value   
#        n.mc.samp,  number of samples from gamma|Y 
#        tau: a set of inverse quantiles corresponding to the edges of the bins of the latent distribution 
#        latent.mixture, input is a list of weights assigned to a uniformly binned latent distribution
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
# Output: a vector of samples from gamma|Y

sample_latent_conditional <- function(y.obs, n.mc.samp, tau, latent.trait.quantiles, cond){
  
  N <- length(cond(0.5)) - 1
  
  # approximate the observed frequencies 
  p.ma <- compute_p_ma(tau, 
                       latent.trait.quantiles,
                       cond)
  
  p.approx <- p.ma[y.obs + 1] # probability of conditional distribution 
  
  # approximate number of times required to sample 
  n.parallel <- round(2*n.mc.samp/p.approx)
  
  # number of bins
  n.bins <- length(tau) - 1
  weights <- c()
  for(i in 1:n.bins){
    weights[i] <- tau[i + 1] - tau[i]
  }
  if(n.parallel >= 500000){
    n.parallel <- 500000
  }
  # outcome 
  out <- c()
  while (length(out) < n.mc.samp){
    
    # sampling from the binned latent distribution
    latent.idx <- sample(1:n.bins, size = n.parallel, replace = TRUE, prob = weights)
    low.bounds <- latent.trait.quantiles[latent.idx]
    high.bounds <- latent.trait.quantiles[latent.idx + 1]
    latent.sample <- runif(n.parallel, min = low.bounds, max = high.bounds)
    
    y.model.samp <- sapply(latent.sample, function(x){
      i = 0:N
      assumption.weights <- cond(x)
      out <- sample(0:N, size = 1, replace = TRUE, prob = assumption.weights)
      return(out)
    })
    
    # keeping the latent variables which generated the observed score 
    keep.latent.idx <- (y.model.samp == y.obs)
    
    out <- c(out, latent.sample[keep.latent.idx])
  }
  # returning only the first exact number of samples 
  out <- out[1:n.mc.samp]
  return(out)
  
}


# Function: sample from the next EM step for updating the latent distribution. 
# Input: y.obs: the observed score value   
#        mu:  regularization parameter
#        tau: a set of inverse quantiles corresponding to the edges of the bins of the latent distribution 
#        latent.mixture: input is a list of weights assigned to a uniformly binned latent distribution
#        cond: function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        n.mc.samp:  number of monte carlo samples from the em update step 
# Output: a vector of samples from gamma|Y

sample_mc_em <- function(p.hat, mu, tau, latent.trait.quantiles, cond, n.mc.samp){
  # tau and quantiles must be of the same length
  # must be ordered 
  
  N <- length(cond(0.5)) - 1
  
  n.tau <- length(tau) - 1
  n.parallel <- n.mc.samp # blocking sampling
  
  weights <- c()
  for(i in 1:n.tau){
    weights[i] <- tau[i + 1] - tau[i]
  }
  # weights based on differences of quantiles 
  
  gamma.em.step <- c()
  
  while (length(gamma.em.step) < n.mc.samp){
    y.hat.samp <- sample(0:N, size = n.parallel, replace = TRUE, prob = p.hat)
    
    frequencies <- n.parallel*compute_edf(y.hat.samp,N)
    # monte carlo sampled from the nonparametric em algorithm 
    gamma.em.step.samp <- lapply(0:N, function(x){
      num.samples <- frequencies[x + 1]
      
      if(num.samples > 0){
        out <- sample_latent_conditional(y.obs = x, 
                                         n.mc.samp = num.samples, 
                                         tau = tau, 
                                         latent.trait.quantiles = latent.trait.quantiles,
                                         cond = cond)
      } else {
        out <- c()
      }
      
      return(out)
    })
    gamma.em.step.samp <- unlist(gamma.em.step.samp)
    gamma.em.step <- c(gamma.em.step, gamma.em.step.samp)
    
  }
  gamma.em.step <- gamma.em.step[1:n.mc.samp]
  
  if(mu > 0){
    n.mixture <- round(n.mc.samp*mu) # round to nearest integer 
    uniform.mixture <- runif(n.mixture)
    out <- c(gamma.em.step,uniform.mixture)
  } else {
    out <- gamma.em.step
  }
  return(out)
}




# Function: computes the loglikelihood value from  
# Input: y.obs: the observed score value   
#        mu:  regularization parameter
#        tau: a set of inverse quantiles corresponding to the edges of the bins of the latent distribution 
#        latent.mixture: input is a list of weights assigned to a uniformly binned latent distribution
#        cond: function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        n.mc.samp:  number of monte carlo samples from the em update step 
# Output: a vector of samples from gamma|Y


compute_loglikelihood_from_latent <- function(p.hat, p.ma, tau, 
                                              latent.trait.quantiles,
                                              mu = 0){
  
  R_bins <- length(tau) - 1
  
  heights <- c()
  for(i in 1:R_bins){
    heights[i] <- (tau[i + 1] - tau[i])/(latent.trait.quantiles[i + 1] - latent.trait.quantiles[i])
  }
  
  widths <- c()
  for(i in 1:R_bins){
    widths[i] <- latent.trait.quantiles[i + 1] - latent.trait.quantiles[i]
  }
  
  if(mu > 0){
    reg.term <- mu*sum(log(heights)*widths)
  } else {
    reg.term <- 0
  }
  
  
  data.term.vec <- p.hat*log(p.ma)
  # defining nan (0*log(0)) values as 0
  data.term.vec[is.nan(data.term.vec)] <- 0 
  
  out <- sum(data.term.vec) + reg.term
  
  return(out)
  
}





# TODO: create a new object which includes the quantiles and the corresponding levels, this will also be the output of this function
fit_nonpar_em <- function(p.hat, cond,  R_bins = 300, mu = 0, n.mc.samp = 1000, verbose = T){
  N <- length(cond(0.5)) - 1
  # looking at a uniform spacing of the inverse quantiles. 
  tau <- seq(0,1,length.out = R_bins)
  # initialize with a uniform distribution
  latent.trait.quantiles.init <- tau
  
  # likelihood change threshold 
  threshold <- 10^(-6)
  
  
  # initialize with a uniform distribution
  current.quantiles <- latent.trait.quantiles.init
  
  n.quantiles <- length(current.quantiles)
  
  # computing the implied distribution on the observed data 
  p.ma <- compute_p_ma(tau = tau,
                       latent.trait.quantiles = current.quantiles,
                       cond = cond,
                       numeric.points = 100)
  
  
  prev.likelihood <- compute_loglikelihood_from_latent(p.hat = p.hat, 
                                                       p.ma = p.ma, 
                                                       tau = tau, 
                                                       latent.trait.quantiles = current.quantiles,
                                                       mu = mu)

  diff.likelihood <- Inf
  em.steps <- 1
  
  while(diff.likelihood >= threshold){
    latent.sample <- sample_mc_em(p.hat = p.hat, 
                                  mu = mu, 
                                  tau = tau, 
                                  latent.trait.quantiles = current.quantiles, 
                                  cond = cond, 
                                  n.mc.samp = n.mc.samp)
    
    # computing updated set of quantiles 
    current.quantiles <- quantile(latent.sample, tau, names = FALSE)
    current.quantiles[1] <- 0
    current.quantiles[n.quantiles] <- 1
    
    p.ma <- compute_p_ma(tau = tau,
                         latent.trait.quantiles = current.quantiles,
                         cond = cond,
                         numeric.points = 100)
    
    new.likelihood <- compute_loglikelihood_from_latent(p.hat = p.hat, 
                                                        p.ma = p.ma, 
                                                        tau = tau, 
                                                        latent.trait.quantiles = current.quantiles,
                                                        mu = mu)
    
    diff.likelihood <- new.likelihood - prev.likelihood
    prev.likelihood <- new.likelihood
    if(verbose){
      print(paste0("EM Steps: ", em.steps, " || Likelihood Change: ", diff.likelihood))
    }
    em.steps <- em.steps + 1
  }
  if(verbose){
    print(paste0("Stochastic EM converged"))
  }
  return(current.quantiles)
}





# Function: compute the A matrix for a binned latent approximation 
# Input: R_bins, the number of uniform latent bins  must an integer  
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        numeric.points,  number of points used in the numeric approximation of A
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value 


compute_A_matrix <- function(R_bins, cond, numeric.points = 100, verbose = F){
  N <- length(cond(0.5)) - 1
  A_matrix <- matrix(data = NA, nrow = N + 1, ncol = R_bins)
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    A.row <- rep(0, N+1)
    for(j in 1:numeric.points){
      y = 0:N
      weights <- cond(design.points[j])
      A.row <- A.row + weights/(sum(weights))
    }
    A.row <- A.row/sum(A.row)
    A_matrix[,i] = A.row
    if(verbose){
      cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
    }
  }
  return(A_matrix)
}



# Function: compute the A tensor using a binned latent approximation 
# Input: R_bins, the number of uniform latent bins  must an integer  
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        numeric.points,  number of points used in the numeric approximation of A
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value 

compute_A_two_obs_tensor <- function(R_bins, cond, numeric.points = 100, verbose = F){
  # 3d Array for faster computation. 
  N <- length(cond(0.5)) - 1
  A_3D <- array(NA, dim = c(N+1,N+1,R_bins))
  
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    
    A.block <- array(0, dim = c(N+1,N+1))
    y1 = 0:N
    y2 = 0:N
    
    for(j in 1:numeric.points){
      
      weights <- outer(cond(design.points[j]), cond(design.points[j]), "*")
      A.block <- A.block + weights/(sum(weights))
    }
    A_3D[,,i] <- A.block/(sum(A.block))
    if(verbose){
      cat(paste0("A Tensor Computed Row: ", i,"/",R_bins), end="\r")
    }
  }
  
  return(A_3D)
}





# Function: Estimate the latent distribution using a convex solver 
# Input: p.hat, a discrete empirical distribution function;  numeric  
#        A.matrix, matrix which maps a vector of uniform bins to the observed data; matrix 
#        mu,  regularization parameter; >= 0
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value 

estimate_mixing_numeric <- function(p.hat, A.matrix, mu){
  
  R.bins <- ncol(A.matrix) # number of bins in the latent distribution
  N <- nrow(A.matrix) - 1 # number of scores 
  theta <- Variable(R.bins, name = "latent discretized distribution") # values of the weight vector 
  obs.dist <- Variable(N + 1, name = "model observed distribution") # values of the observed distribution
  
  data.obj <-  t(p.hat) %*% log(obs.dist)
  pen.obj <- (mu/R.bins)*t(rep(1, R.bins)) %*% log(theta) 
  
  constraints <- list(
    obs.dist == A.matrix %*% theta,
    sum(theta) <= 1,
    mu/(R.bins*(1 + mu)) <= theta
  )
  
  obj.arg <- data.obj + pen.obj
  obj <- Maximize(obj.arg)
  p <- Problem(obj, constraints)
  
  
  #value(theta) <- rep(1/R.bins, R.bins) # initial guess of a uniform distribution
  result <- solve(p, solver = "MOSEK")
  #result <- solve(p, verbose = TRUE)
  
  p.m <- result$getValue(theta)
  p.ma <- result$getValue(obs.dist)
  out.list <- list("latent" = p.m, "observed" = p.ma)
  return(out.list)
  
}


##   Feasibility tests ---------

# Function: computes the first order feasibility test 
# Input: latent.mixture, a vector of weights for the latent distribution; numeric            
#        A.matrix, matrix which maps a vector of uniform bins to the observed data; matrix 
#        p.hat, a discrete empirical distribution function;  must be an vector of numeric values  
# Output: p.feasibility: p value for the feasibility test, is at minimum e-8 due to numerical rounding issues 

test_feasibility_first_order <- function(latent.mixture, A.matrix, p.hat, sample.size){
  # model implied distribution on Y 
  p.ma <- A.matrix %*% latent.mixture
  p.ma <- as.numeric(p.ma)
  # likelihood ratio statistic
  lr <- 2*sample.size*kl_divergence(p.hat, p.ma)
  k <- length(p.hat)
  
  # threshold for a very low p.value. 
  # numerical errors can occur if the likelihood ratio statistic is too large
  # for values which would correspond to p-values below 10^(-8) we simply use 10^(-8) 
  lr.thresh <- criticalValue(k,sample.size, p = 10^(-8))
  if(lr >= lr.thresh){
    p.feasibility <- 10^(-8)
  } else {
    p.feasibility <- tailProbBound(x = lr, k = k, n = sample.size)
    if(p.feasibility > 1){
      p.feasibility <- 1
    }
    # handling rounding errors  
    if(p.feasibility < 0){
      p.feasibility <- 0
    }
  }
  return(p.feasibility)
}


# Function: computes the second order feasibility test 
# Input: latent.mixture, a vector of weights for the latent distribution; numeric            
#        A.two.sample.tensor, tensor which maps a vector of uniform bins to the bivariate observed data distribution; 3D tensor 
#        p.hat, a discrete empirical distribution function;  must be an vector of numeric values  
# Output: p.feasibility: p value for the feasibility test, is at minimum e-8 due to numerical rounding issues 

test_feasibility_second_order <- function(latent.mixture, A.two.sample.tensor, p.hat, sample.size){
  # model implied distribution on Y 
  R_bins <- length(latent.mixture)
  p.ma.list <- lapply(1:R_bins, function(z){
    out <- A.two.sample.tensor[,,z]*latent.mixture[z]
    return(out)
    })
  p.ma <- matrix(data = 0, nrow = nrow(A.two.sample.tensor[,,1]),
                 ncol = ncol(A.two.sample.tensor[,,1]))
  
  for(i in 1:R_bins){
    p.ma <- p.ma + p.ma.list[[i]]
  }
  
  # likelihood ratio statistic
  lr <- 2*sample.size*kl_divergence(p.hat, p.ma)
  
  k <- nrow(p.hat)*ncol(p.hat)
  # threshold for a very low p.value. 
  # numerical errors can occur if the likelihood ratio statistic is too large
  # for values which would correspond to p-values below 10^(-8) we simply use 10^(-8) 
  lr.thresh <- criticalValue(k,sample.size, p = 10^(-8))
  if(lr >= lr.thresh){
    p.feasibility <- 10^(-8)
  } else {
    p.feasibility <- tailProbBound(x = lr, k = k, n = sample.size)
    if(p.feasibility > 1){
      p.feasibility <- 1
    }
    # handling rounding errors  
    if(p.feasibility < 0){
      p.feasibility <- 0
    }
  }
  return(p.feasibility)
}



##   Model Selection Functions ----------------------------------------

# Function: Numerically approximate the distribution of second test score, given the first
# Input: y.obs, the observed score;  must be an integer  
#        latent.mixture, weights assigning to the discretized latent distribution; vector summing to 1
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> P{0, ..., N}
# Output: probability : weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value 



second_score_conditional <- function(y.obs, latent.mixture, cond){
  
  tau <- inv_quantiles_from_weights(latent.mixture)
  
  # corresponding quantiles to the latent mixture
  latent.quantiles <- seq(0,1, length.out = length(tau))
  
  R_bins <- length(latent.mixture)
  latent.points <- sapply(1:R_bins, function(z){
    gam <- (latent.quantiles[z] + latent.quantiles[z+1])/2
    return(gam)
  })
  
  
  joint.dist.grid  <- sapply(1:R_bins, function(gam.idx){
    gam <- latent.points[gam.idx]
    weight <- latent.mixture[gam.idx]
    
    out.y1y2.joint <- cond(gam) %o% cond(gam)
    prob.row <- out.y1y2.joint[,y.obs + 1]
    out <- prob.row*weight
    return(out)
  })
  
  cond.prob <- rowSums(joint.dist.grid)
  cond.prob <- cond.prob/sum(cond.prob)
  return(cond.prob)
}



# TODO: Edit this lil guy here 
# Function: Estimate the latent distribution using a convex solver 
# Input: p.hat, a discrete empirical distribution function;  numeric  
#        A.matrix, matrix which maps a vector of uniform bins to the observed data; matrix 
#        mu,  regularization parameter; >= 0
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value 


compute_two_obs_loglikelihood <- function(y.paired, A.two.sample.tensor, latent.mixture.list){
  M <- nrow(y.paired)
  log.likelihood <- 0
  for(m in 1:M){
    pair.prob <- t(A.two.sample.tensor[y.paired[m,1] +1,y.paired[m,2] + 1,]) %*% latent.mixture.list[[m]]
    log.likelihood <- log.likelihood + as.numeric(log(pair.prob))
  }
  return(log.likelihood)
}


# TODO: document this 
compute_joint_dist_beta <- function(beta1.model.y, beta2.model.y,
                                    beta1.model.z, beta2.model.z,
                                    cond.y, cond.z, grid.size = 10000){
  
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  
  latent.grid <- seq(0,1, length.out = grid.size)
  
  joint.dist <- lapply(latent.grid, function(x){
    latent.y <- qbeta(x, shape1 = beta1.model.y, shape2 = beta2.model.y)
    latent.z <- qbeta(x, shape1 = beta1.model.z, shape2 = beta2.model.z)
    iy = 0:Ny
    iz = 0:Nz
    assumption.weights.y <- cond.y(latent.y)
    assumption.weights.z <- cond.z(latent.z)
    
    prob.y <- assumption.weights.y/sum(assumption.weights.y)
    prob.z <- assumption.weights.z/sum(assumption.weights.z)
    
    p.yz <- prob.y %*% t(prob.z)
    out <- p.yz
    
  })
  
  p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
  for(i in 1:length(joint.dist)){
    p.yz <- p.yz + joint.dist[[i]]
  }
  p.yz <- p.yz/grid.size
  return(p.yz)
}






##### TO DO: EDIT THIS FUNCTION
# Function: Sample from the model implied distribution on a subsequent score assuming the latent variable remains constant 
# Input: pair.obs, an observed data set of pairs of distributions for which
#                  the latent variable is constant;  must be an vector of integers  
#        latent.mixture, list of latent mixtures ; input is a list of weights assigned to a uniformly binned latent  
#        
# Output: a list of true.diff: The true difference between subsequent scores 
#                    model.sample.diff: The difference between the model sampled second score and the true first score


intrinsic_variability_sample <- function(pair.obs, n.samp, latent.mixture, cond){
  N <- length(cond(0.5)) - 1
  # sampling from the model implied version of the subsequent conditional 
  y.model.sample <- sample(0:N, size = n.samp, replace = T, prob = second_score_conditional(pair.obs[1], latent.mixture, cond))
  
  d.true <- pair.obs[2] - pair.obs[1]
  d.sim <- y.model.sample - pair.obs[1]
  
  out.list <- list("true.diff" = d.true, "model.sample.diff" = d.sim)
  return(out.list)
}





# Function: compute the intrinsic variability from a data sample of pair.obs
# Input: pair.obs, an observed data set of pairs of distributions for which
#                  the latent variable is constant;  must be an vector of integers  
#        latent.mix.list, list of latent mixtures ; input is a list of weights assigned to a uniformly binned latent  
#        model.observed.list,  list of corresponding distributions on the observed data
# Output: latent: weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value 
# TODO: Delete model.observed.list

intrinsic_variability <- function(pair.obs, latent.mix.list, n.samp, cond){
  # samples from the true q0 distribution 
  e.true <- pair.obs[,2] - pair.obs[,1]
  n <- nrow(pair.obs)
  
  e.sim.list <- lapply(1:n, function(x){
    latent.mix <- latent.mix.list[[x]]
    intrinsic.samp <- intrinsic_variability_sample(pair.obs = as.numeric(pair.obs[x,]), 
                                                   n.samp = n.samp, latent.mixture = latent.mix, 
                                                   cond = cond)
    return(intrinsic.samp$model.sample.diff)
  })
  e.sim <- unlist(e.sim.list)
  # intrinsic variability and the norm 
  int.norm <- tv_norm(e.true,e.sim)
  return(int.norm)
}









##   Conversion Functions ----------------------------------------





#### Score Conversion Functions 

### TO DO: EDIT THIS FUNCTION description 
# Function: Numerically approximate the distribution of second scores, given the first
# Input: y.obs, the observed score;  must be an integer  
#        latent.mixture, weights assigning to the discretized latent distribution; vector summing to 1
#        cond, function defining the conditional distribution of Y|gamma; function [0,1] -> P{0, ..., N}
# Output: probability : weights in each bin in the latent distribution
#         observed: list of probabilities assigned to each test score value 


compute_conversion_prob <- function(y, z, latent.mixture.y, latent.mixture.z, 
                                    cond.y, cond.z, grid.size = 10000){
  
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  
  tau.y <- inv_quantiles_from_weights(latent.mixture.y)
  latent.quantiles.y <- seq(0,1, length.out = length(tau.y))
  
  tau.z <- inv_quantiles_from_weights(latent.mixture.z)
  latent.quantiles.z <- seq(0,1, length.out = length(tau.z))
  
  latent.grid <- seq(0,1, length.out = grid.size)
  
  joint.dist <- lapply(latent.grid, function(x){
    latent.y <- qf_L(x, latent.quantiles.y, tau.y)
    latent.z <- qf_L(x, latent.quantiles.z, tau.z)
    prob.y <- cond.y(latent.y)
    prob.z <- cond.z(latent.z)
    out <- prob.y %*% t(prob.z)
  })
  
  p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
  for(i in 1:length(joint.dist)){
    p.yz <- p.yz + joint.dist[[i]]
  }
  p.yz <- p.yz/grid.size
  
  pz.given.y <- exp(log(p.yz[y + 1,] ) - log(sum(p.yz[y + 1,])))
  
  cond.prob <- pz.given.y[z + 1]
  return(cond.prob)
}




# TODO: name and document 
# replace _metric with _ce 
convert_score_metric <- function(test.pairs, latent.mix.list.y, latent.mix.list.z, 
                                 cond.y, cond.z, joint.prob, grid.size = 1000){
  
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  
  # weighting for whether we have access to the joint probabilities 
  
  if(missing(joint.prob)){
    joint.prob <- matrix(data = 1, nrow = Ny + 1, ncol = Nz + 1)
  }
  
  z.true <- test.pairs[,2]
  z.sim <- c()
  
  # TODO: problem indices
  idx <- 1:nrow(test.pairs)
  res <- sapply(idx, function(x){
    latent.mix.y <- latent.mix.list.y[[x]]
    latent.mix.z <- latent.mix.list.z[[x]]
    
    
    
    conv.prob <- compute_conversion_prob(y = test.pairs[x,1], 
                                         z = test.pairs[x,2], 
                                         latent.mixture.y = latent.mix.y, 
                                         latent.mixture.z = latent.mix.z, 
                                         cond.y = cond.y, 
                                         cond.z = cond.z, 
                                         grid.size = grid.size)
    
    
    joint.weight <- joint.prob[test.pairs[x,1] + 1, test.pairs[x,2] + 1]
    result.part <- -joint.weight*log(conv.prob)
    
    
    return(result.part)
  })
  
  out <- sum(res)
  
  return(out)
}




##   Parametric model functions  ---------------------------------------------

# TODO: Document and delete the rest
logitnorm_model_fit <- function(outcome, X, cond, numeric.points = 1000){
  N = length(length(cond(0.5))) - 1
  renorm.r.function.denom <- rep(0,N + 1)
  numeric.grid <- seq(0,1,length.out = numeric.points)
  for(x in numeric.grid){
    renorm.r.function.denom <- renorm.r.function.denom + cond(x)
  }
  mean.seq <- rep(0, N + 1)
  mom2.seq <- rep(0, N + 1) # 2nd moment sequence
  trim.grid <- numeric.grid[2:(length(numeric.grid) - 1)] # removes problematic 0 and 1 
  for(x in trim.grid){
    mean.seq <- mean.seq + cond(x)*logit(x)/renorm.r.function.denom
  }
  for(x in trim.grid){
    mom2.seq <- mom2.seq + cond(x)*logit(x)^2/renorm.r.function.denom
  }
  
  y <- mean.seq[outcome + 1]
  
  beta.coef <- solve(t(as.matrix(X)) %*% as.matrix(X)) %*% (t(as.matrix(X)) %*% y)
  
  n <- length(y)
  mus <-  as.matrix(X) %*% beta.coef
  sig2 <- mean(mom2.seq[outcome + 1] -2*mus*y +mus^2)
  sigma <- sqrt(sig2)
  
  out.params <- list("beta" = beta.coef, "sigma" = sigma)
  return(out.params)
}

compute_conversion_parametric <- function(y, z, X.row, params.y, params.z, cond.y, cond.z, grid.size = 1000){
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  
  mu.y <- as.numeric(X.row) %*% params.y$beta
  mu.z <- as.numeric(X.row) %*% params.z$beta
  sig.y <- params.y$sigma
  sig.z <- params.z$sigma
  
  latent.grid <- seq(0,1, length.out = grid.size)
  #trim.grid <- latent.grid[2:(length(latent.grid) - 1)] # removes the troublesome end points 
  
  # here we compute the model implied joint probability using a grid 
  # of points and the inverse quantile function of the normal. 
  joint.dist <- lapply(latent.grid, function(x){
    latent.y <- logistic(qnorm(x, mean = mu.y, sd = sig.y))
    latent.z <- logistic(qnorm(x, mean = mu.z, sd = sig.z))
    
    assumption.weights.y <- cond.y(latent.y)
    assumption.weights.z <- cond.z(latent.z)
    
    prob.y <- assumption.weights.y/sum(assumption.weights.y)
    prob.z <- assumption.weights.z/sum(assumption.weights.z)
    out <- prob.y %*% t(prob.z)
  })
  
  p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
  for(i in 1:length(joint.dist)){
    p.yz <- p.yz + joint.dist[[i]]
  }
  p.yz <- p.yz/grid.size
  
  pz.given.y <- p.yz[y + 1,]/sum(p.yz[y + 1,])
  
  cond.prob <- pz.given.y[z + 1]
  return(cond.prob)
}



compute_conversion_parametric_ce <- function(test.pairs, X, params.y, params.z, 
                                             cond.y, cond.z, grid.size = 1000){
  
  Ny <- length(cond.y(0.5)) - 1
  Nz <- length(cond.z(0.5)) - 1
  
  # weighting for whether we have access to the joint probabilities 
  
  
  z.true <- test.pairs[,2]
  z.sim <- c()
  
  
  idx <- 1:nrow(test.pairs)
  res <- sapply(idx, function(x){
    conv.prob <- compute_conversion_parametric(test.pairs[x,1], test.pairs[x,2],  
                                               X.row = X[x,], 
                                               params.y = params.y, params.z = params.z, 
                                               cond.y = cond.y, 
                                               cond.z = cond.z, 
                                               grid.size = grid.size)
    
    
    
    result.part <- -log(conv.prob)
    
    
    return(result.part)
  })
  
  out <- sum(res)
  
  return(out)
}







# TODO: rename and clarify what this is for (parametric)
renormalized_r_function <- function(N, ker, h, is.binomial = F, grid.size = 100000){
  norm_consts <- rep(0,N+1)
  eps = 10^(-6)
  y.vals <- 0:N
  gam.grid <- seq(0,1, length.out = grid.size)
  trim.grid <- seq(eps,1 - eps, length.out = grid.size)
  mean.seq <- rep(NA, N + 1)
  mom2.seq <- rep(NA, N + 1)
  if(is.binomial){
    if(missing(ker) | missing(h)){
      ker <- NA
      h <- NA
    }
    for(k in y.vals){
      norm_consts[k + 1] = mean(dbinom(k, size = N, prob = gam.grid))
      mean.seq[k + 1] = mean(dbinom(k, size = N, prob = trim.grid)*logit(trim.grid))/norm_consts[k + 1]
      mom2.seq[k + 1] = mean(dbinom(k, size = N, prob = trim.grid)*logit(trim.grid)^2)/norm_consts[k + 1]
    }
    
  }
  else{
    for(k in y.vals){
      norm_consts[k + 1] = mean(ker((k - gam.grid*N)/h))
      mean.seq[k + 1] = mean(ker((k - gam.grid*N)/h)*logit(trim.grid))/norm_consts[k + 1]
      mom2.seq[k + 1] = mean(ker((k - gam.grid*N)/h)*logit(trim.grid)^2)/norm_consts[k + 1]
    }
  }
  out <- list("latent.logit.mean" = mean.seq, "latent.logit.mom2" = mom2.seq)
  return(out)
}


# TODO: rename, make sure its obvious about the latent version 
fit_logitnorm_model <- function(k, X, latent.variational.moments){
  # Ensure that we have a matrix/ dataframe and a density object. 
  # should also ensure that it had been trained
  
  
  y <- latent.variational.moments$latent.logit.mean[k + 1]
  
  beta_coef <-solve(t(as.matrix(X)) %*% as.matrix(X)) %*% (t(as.matrix(X)) %*% y)
  # Could probably paralellize the variance computation 
  n <- length(y)
  mus <-  as.matrix(X) %*% beta_coef
  sig2 <- mean(latent.variational.moments$latent.logit.mom2[k + 1] -2*mus*latent.variational.moments$latent.logit.mean[k + 1] +mus^2)
  sigma <- sqrt(sig2)
  
  data_out <- list(beta = beta_coef, sigma = sigma)
  return(data_out)
  
}



### TODO: Change this and stick to functional programming 

# TODO: rename, to psi function, parametric conversion 
conv_fx <- setRefClass("Parametric Conversion Function",
                       fields = list(y_params = "list", z_params = "list"),
                       methods = list(
                         convert_yz = function(gam,x){
                           zeta <- logistic(z_params$sigma/y_params$sigma *logit(gam) + 
                                              sum((z_params$beta - (z_params$sigma/y_params$sigma)*y_params$beta)*x))
                           return(zeta)
                         },
                         convert_zy = function(zeta,x) {
                           gam <- logistic(y_params$sigma/z_params$sigma *logit(zeta) + 
                                             sum((y_params$beta - (y_params$sigma/z_params$sigma)*z_params$beta)*x))
                           return(gam)
                         }
                       )
)

# TODO: rename, make sure it is obvious that this is an estimator of the conditional z|y 
# TODO: Parametric Feels a bit like filler
conversion.parametric <- function(y, z, x, params.y, params.z, 
                                  Ny, Nz, ker.y, ker.z, hy, hz, 
                                  binomial.y = F, binomial.z = F,  R_bins = 1000){
  
  if((missing(ker.y) | missing(hy) ) & ! binomial.y){
    stop("Missing Meaurement model on Y")
  }
  
  if((missing(ker.z) | missing(hz) ) & ! binomial.z){
    stop("Missing Meaurement model on Z")
  }
  
  if((missing(ker.y) | missing(hy))){
    ker.y <- NA
    hy <- NA
  }
  
  if((missing(ker.z) | missing(hz))){
    ker.z <- NA
    hz <- NA
  }
  
  
  
  
  ######################################################
  ######################################################
  
  if(!binomial.y){
    p.ay <- function(y,gam){
      i <- 0:Ny
      out <- ker.y((y - Ny*gam)/hy)/sum(ker.y((i - Ny*gam)/hy))
      return(out)
    }
  }
  
  
  if(!binomial.z){
    p.az <- function(z,zeta){
      j <- 0:Nz
      out <- ker.z((z - Nz*zeta)/hz)/sum(ker.z((j - Nz*zeta)/hz))
      return(out)
    }
  }
  
  if(binomial.y){
    p.ay <- function(y,gam){
      out <- dbinom(x = y, size = Ny, prob = gam)
      return(out)
    }
  }
  
  
  if(binomial.z){
    p.az <- function(z,zeta){
      out <- dbinom(x = z, size = Nz, prob = zeta)
      return(out)
    }
  }
  
  
  mu.y <- t(params.y$beta) %*% as.numeric(x) 
  mu.z <- t(params.z$beta) %*% as.numeric(x) 
  
  sigma.y <- params.y$sigma
  sigma.z <- params.z$sigma
  
  conversion <- conv_fx(y_params = params.y, z_params = params.z)
  
  ######
  # plot(latent.grid, conversion$convert_yz(latent.grid, x))
  #####
  
  eps = 10^(-8)
  latent.grid <-seq(eps, 1 - eps, length.out = R_bins) 
  joint.dist.grid  <- sapply(1:R_bins, function(gam.idx){
    gam <- latent.grid[gam.idx]
    weight <- 1/(gam*(1 - gam))*dnorm(logit(gam), mean = mu.y, sd = sigma.y)
    zeta.point <- conversion$convert_yz(gam, x)
    
    out.y <- p.ay(y,gam)
    out.z <- p.az(z,zeta.point)
    out <- out.y*out.z*weight
    return(out)
  })
  
  # density function of a logit-normal 
  p.m <- 1/(latent.grid*(1 - latent.grid))*dnorm(logit(latent.grid), mean = mu.y, sd = sigma.y)
  denom <- mean(p.ay(y, latent.grid)*p.m)
  cond.out <- mean(joint.dist.grid)/(denom)
  return(cond.out)
}


# TODO: should not also need this
par.cond.sample <- function(y,N,mu,sig, p.ay, n.samp){
  n_parallel <- 100000
  gam.out <- c()
  while(length(gam.out) < n.samp ){
    norm.samp <- rnorm(n_parallel, mean = mu, sd = sig)
    logit.norm.samp <- logistic(norm.samp)
    y.samp <- sapply(logit.norm.samp, function(z){
      out <- sample(0:N, size = 1, prob = p.ay(0:N, z), replace = T)
      return(out)
    })
    keep.idx <- y.samp == y
    gam.out <- c(gam.out, logit.norm.samp[keep.idx])
  }
  return(gam.out)
}


# TODO: rename and change the . to _ 
naive.conversion.prob <- function(y,muy,sdy,muz,sdz, Nz){
  z.pred <- muz + (sdz/sdy)*(y - muy)
  
  
  
  z.scale <- 0:Nz
  
  out.prob <- sapply(z.scale, function(z){
    z.max <- z + 1/2
    z.min <- z - 1/2
    if(z.max > Nz){
      z.max = Inf
    }
    if(z.min < 0){
      z.min = -Inf
    }
    
    out <- pnorm(z.max, mean = z.pred, sd = sdz) - pnorm(z.min, mean = z.pred, sd = sdz)
    return(out)
  })
  
  return(out.prob)
}







##   Functions for simulations ---------------------------------------------



# Function: simulated a latent beta distribution
# Input: n.ind, number of individual;  must be an integer 
#        n.obs.per.ind, number of repeated observations per individual; integer 
#        (beta1.y, beta2.y), parameters for the beta distribution in the model; > 0
#        cond.y, function defining the conditional distribution of Y|gamma; function [0,1] -> {0, ..., N}
#        pair.obs, logical whether sampling from 
# Output: data set of sampled according to the generated model 
simulate_beta <- function(n.ind, n.obs.per.ind, beta1.y, beta2.y, cond.y, beta1.z, beta2.z, cond.z, pair.obs = T){
  
  # simulating from the model 
  omega <- runif(n.ind, 0, 1)
  latent.samp.y <- qbeta(omega, beta1.y, beta2.y)
  y.columns <- sprintf("y%s",seq(1:n.obs.per.ind))
  Ny <- length(cond.y(0.5)) - 1 
  if(pair.obs) {
    if(missing(beta1.z) | missing(beta2.z) | missing(cond.z)){
      stop("Must include model on Z branch")
    }
    
    latent.samp.z <- qbeta(omega, beta1.z, beta2.z)
    z.columns <- sprintf("z%s",seq(1:n.obs.per.ind))
    Nz <- length(cond.z(0.5)) - 1 
    
    dataset <- data.frame(matrix(NA, nrow = n.ind, ncol = 2*n.obs.per.ind))
    colnames(dataset) <- c(y.columns,z.columns)
    # Sample jointly from a fixed latent variable 
    for(j in 1:n.obs.per.ind){
      y.samp <- sapply(latent.samp.y, function(s){
        out <- sample(0:Ny,1, replace = T, prob = cond.y(s))
        return(out)
      })
      dataset[,j] <- y.samp
      
      z.samp <- sapply(latent.samp.z, function(s){
        out <- sample(0:Nz,1, replace = T, prob = cond.z(s))
        return(out)
      })
      dataset[,j + n.obs.per.ind] <- z.samp
    }
  } else{ # Case where only y is sampled 
    dataset <- data.frame(matrix(NA, nrow = n.ind, ncol = n.obs.per.ind))
    colnames(dataset) <- c(y.columns)
    # Sample jointly from a fixed latent variable 
    for(j in 1:n.obs.per.ind){
      y.samp <- sapply(latent.samp.y, function(s){
        out <- sample(0:Ny,1, replace = T, prob = cond.y(s))
        return(out)
      })
      dataset[,j] <- y.samp
    }
  }
  return(dataset)
}





#TODO: Document 
summarize_simulations <- function(results){
  
  array.shape <- dim(results)
  n.samples <- array.shape[1]
  n.dim <- length(array.shape)
   
  
  if(n.dim > 2 ){
    mean.array <- array(0, dim = array.shape[2:length(array.shape)])
    for(i in 1:n.samples){
      sim.array <- extract(results, "1"=i)
      sim.array <- wrap(sim.array, map=as.list(c(NA,3:n.dim)))
      mean.array <- mean.array + sim.array
    }
    mean.array <- mean.array/n.samples
    
    var.array <- array(0, dim = array.shape[2:length(array.shape)]) 
    
    for(i in 1:n.samples){
      sim.array <- extract(results, "1"=i)
      sim.array <- wrap(sim.array, map=as.list(c(NA,3:n.dim)))
      var.array <- var.array + (mean.array - sim.array)^2
    }
    var.array <- var.array/(n.samples - 1)
    
    # standard error of the mean estimate (assuming LLN)
    s.error.array <- sqrt(var.array/n.samples)
  } else if (n.dim == 2){
    mean.array <- array(0, dim = array.shape[2:length(array.shape)])
    for(i in 1:n.samples){
      sim.array <- results[i,]
      mean.array <- mean.array + sim.array
    }
    mean.array <- mean.array/n.samples
    
    var.array <- array(0, dim = array.shape[2:length(array.shape)]) 
    
    for(i in 1:n.samples){
      sim.array <- results[i,]
      var.array <- var.array + (mean.array - sim.array)^2
    }
    var.array <- var.array/(n.samples - 1)
    
    # standard error of the mean estimate (assuming LLN)
    s.error.array <- sqrt(var.array/n.samples)
    
  } else {
    stop("array must be have at least 2 dimensions")
  }

  out <- list("mean_results" = mean.array, 
              "se_results" = s.error.array)
  
  return(out)
  
}



##   Data Analysis functions  ------------------------------------------------------


# Function: groups the sex and education variables into 4 particular categories (used for cleaning the real data set)
# Input: data set with sex and educ as variables 
# Output: new data set replacing those two with a corresponding group 
#       1: female, low education 
#       2: male, low education 
#       3: female, high education 
#       4: male, high education 

categorize <- function(data){
  tmp <- data %>% mutate(group = ifelse(((sex == 1) & (educ <= 16)), 1, 
                                        ifelse(((sex == 2) & (educ <= 16)), 2, 
                                               ifelse(((sex == 1) & (educ > 16)), 3, 4))))
  tmp <- tmp[, !names(tmp) %in% c('sex','educ')]  
  tmp$group <- factor(tmp$group)
  return(tmp)
}


### TODO: Document 


naive_conversion_prob <- function(y,muy,sdy,muz,sdz, Nz, n.samp = 100000){
  z.pred <- muz + (sdz/sdy)*(y - muy)
  eps <- rnorm(n.samp, mean = 0, sd = sdz)
  
  z.samp <- z.pred + eps
  
  z.scale <- 0:Nz
  
  closest.idx <- sapply(z.samp, function(z){
    dists <- abs(z.scale - z)
    out <- which.min(dists)
    return(out)
  })
  
  rounded.samp <- z.scale[closest.idx]
  out <- rep(0,Nz + 1)
  for(k in 0:Nz){
    out[k + 1] <- mean(rounded.samp == k)
  }
  
  return(out)
}



##   TRASHBIN ---------------------------------------------
# TO DO: delete after using simplified one 
cond.dist.est <- function(x.design, train.data, outcome, N, age.window = 3){
  out.p.hat <- matrix(data = NA, nrow = nrow(x.design), ncol = N + 1)
  for(i in 1:nrow(x.design)){
    x <- x.design[i,]
    idx1 <- x$group == train.data$group 
    idx2 <- abs(x$age - train.data$age) <= age.window
    idx <- idx1 & idx2
    sub.dat <- train.data[idx,]
    obs.set <- 0:N
    summary.dat <- sub.dat %>%
      group_by_(outcome) %>%
      summarise(Count = n()) 
    summary.dat <- as.data.frame(summary.dat)
    
    n <- sum(summary.dat$Count)
    
    missing.obs <- obs.set[!obs.set %in% summary.dat[1:nrow(summary.dat),outcome]]
    missing.block <- data.frame(y = missing.obs, Count = rep(0,length(missing.obs)))
    colnames(missing.block)[1] <- outcome
    summary.dat <- rbind(summary.dat, missing.block)
    
    summary.dat <- summary.dat %>% arrange_(outcome)
    
    p.hat <- summary.dat$Count/sum(summary.dat$Count)
    out.p.hat[i,] <- p.hat
  }
  
  return(out.p.hat)
}


