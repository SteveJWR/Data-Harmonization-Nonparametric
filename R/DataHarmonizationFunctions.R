### 


## TODO: Create categories of functions
# 

logit <- function(x){
  if ({
    any(x < 0)
  } || {
    any(x > 1)
  }) 
    stop("x must be in [0,1].")
  out = log(x/(1-x))
  return(out)
}

logistic <- function(x){
  return(1/(1 + exp(-x)))
}


# TO DO: Create a dictionary for the conditional distributions 
categorize <- function(data){
  tmp <- data %>% mutate(group = ifelse(((sex == 1) & (educ <= 16)), 1, 
                                        ifelse(((sex == 2) & (educ <= 16)), 2, 
                                               ifelse(((sex == 1) & (educ > 16)), 3, 4))))
  tmp <- tmp[, !names(tmp) %in% c('sex','educ')]  
  tmp$group <- factor(tmp$group)
  return(tmp)
}


# TO DO: See if it is possible to run this just importing this from Richard's package

tailProbBound <- function (x, k, n, verbose = FALSE){
  stopifnot(`k and n must be of the same length` = length(k) == 
              length(n))
  
  
  
  if (length(k) == 1) {
    f <- function(lambda) {
      lambda <- pmin(pmax(lambda, 0), 1)
      -lambda * x/2 + log(mgfBound(k, n, lambda))
    }
  }
  else {
    f <- function(lambda) {
      lambda <- pmin(pmax(lambda, 0), 1)
      result <- -lambda * x/2
      for (i in 1:length(k)) {
        result <- result + log(mgfBound(k[i], n[i], lambda))
      }
      return(result)
    }
  }
  if (length(n) == 1) {
    lambda.init <- max(min(1 - (k - 1)/(x/2) + k/(k - 1) * 
                             (x/2 - k + 1)/n, 1), 0)
    lambda.init <- c(lambda.init, seq(0, 1, length.out = 3))
  }
  else {
    lambda.init <- seq(0, 1, length.out = 4)
  }
  
  # additional line for divergences above .5 to infinity 
  lambda.init <- lambda.init[lambda.init < .2]
  
  sols <- plyr::laply(lambda.init, function(.x) {
    stats::optim(.x, f, method = "L-BFGS-B", lower = 0 , upper = 1)
  })
  m <- which.min(sols[, 2])
  sol.par <- unlist(sols[, 1])[m]
  sol.val <- unlist(sols[, 2])[m]
  if (verbose) {
    lambda.vec <- seq(0, 1, length.out = 100)
    plot(lambda.vec, f(lambda.vec), type = "l", xlab = "lambda", 
         ylab = "bound")
    graphics::abline(v = sol.par, col = "red")
  }
  return(exp(sol.val))
}

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



measurement_kernel_sample <- function(size, lat, N, ker, h){
  i = 0:N
  probs <- ker((i - N*lat)/h)
  probs <- probs/(sum(probs))
  out <- sample(0:N, size = size, replace = TRUE, prob = probs)
  return(out)
}


# TODO: rename all sample_latent_conditional
conditional.samp <- function(y.cond, n.samp,  p.hat, tau, latent.trait.quantiles, N, ker, h){
  
  p.approx <- p.hat[y.cond + 1] # probability of conditional distribution 
  
  n_parallel <- round(2*n.samp/p.approx)
  
  n_tau <- length(tau) - 1
  weights <- c()
  for(i in 1:n_tau){
    weights[i] <- tau[i + 1] - tau[i]
  }
  if(n_parallel >= 50000000){
    n_parallel <- 50000000
  }
  # outcome 
  out <- c()
  while (length(out) < n.samp){
    latent.idx <- sample(1:n_tau, size = n_parallel, replace = TRUE, prob = weights)
    low_bounds <- latent.trait.quantiles[latent.idx]
    high_bounds <- latent.trait.quantiles[latent.idx + 1]
    latent.sample <- runif(n_parallel, min = low_bounds, max = high_bounds)
    
    y.model.samp <- lapply(latent.sample, function(x){
      i = 0:N
      assumption.weights <- ker((i - N*x)/h)
      out <- sample(0:N, size = 1, replace = TRUE, prob = assumption.weights)
      return(out)
    })
    y.model.samp <- unlist(y.model.samp)
    
    keep.latent.idx <- (y.model.samp == y.cond)
    
    out <- c(out, latent.sample[keep.latent.idx])
    #print(paste0(length(out), " of ", n.samp))
  }
  out <- out[1:n.samp]
  return(out)
  
}

# TODO: create a version which samples from binomial or the kernel method
conditional.binom.samp <- function(y.cond, n.samp,  p.hat, tau, latent.trait.quantiles, N){
  
  p.approx <- p.hat[y.cond + 1] # probability of conditional distribution 
  
  n_parallel <- round(2*n.samp/p.approx)
  
  n_tau <- length(tau) - 1
  weights <- c()
  for(i in 1:n_tau){
    weights[i] <- tau[i + 1] - tau[i]
  }
  if(n_parallel >= 50000000){
    n_parallel <- 50000000
  }
  # outcome 
  out <- c()
  while (length(out) < n.samp){
    latent.idx <- sample(1:n_tau, size = n_parallel, replace = TRUE, prob = weights)
    low_bounds <- latent.trait.quantiles[latent.idx]
    high_bounds <- latent.trait.quantiles[latent.idx + 1]
    latent.sample <- runif(n_parallel, min = low_bounds, max = high_bounds)
    
    y.model.samp <- sapply(latent.sample, function(x){
      i = 0:N
      out <- rbinom(1, size = N, prob = x)
      return(out)
    })
    
    keep.latent.idx <- (y.model.samp == y.cond)
    
    out <- c(out, latent.sample[keep.latent.idx])
    #print(paste0(length(out), " of ", n.samp))
  }
  out <- out[1:n.samp]
  return(out)
  
}


# TODO: create a version which samples from binomial or the kernel method
# TODO: rename and replace . with _ 
mc.em.samp <- function(p.hat, tau, latent.trait.quantiles, N, ker, J, h){
  # tau and quantiles must be of the same length
  # must be ordered 
  
  
  n_tau <- length(tau) - 1
  n_parallel <- J # blocking sampling
  
  weights <- c()
  for(i in 1:n_tau){
    weights[i] <- tau[i + 1] - tau[i]
  }
  # weights based on differences of quantiles 
  
  out <- c()
  
  while (length(out) < J){
    # latent.idx <- sample(1:n_tau, size = n_parallel, replace = TRUE, prob = weights)
    # low_bounds <- latent.trait.quantiles[latent.idx]
    # high_bounds <- latent.trait.quantiles[latent.idx + 1]
    # latent.sample <- runif(n_parallel, min = low_bounds, max = high_bounds)
    # 
    y.hat.samp <- sample(0:N, size = n_parallel, replace = TRUE, prob = p.hat)
    
    gamma.step.samp <- lapply(y.hat.samp, function(x){
      
      out.tmp <- conditional.samp(y.cond = x, n.samp = 1,  p.hat = p.hat, tau = tau, 
                                  latent.trait.quantiles = latent.trait.quantiles, 
                                  N = N, ker = ker, h = h)
      return(out.tmp)
    })
    gamma.step.samp <- unlist(gamma.step.samp)
    out <- c(out, gamma.step.samp)
    #print(paste0(length(out), " of ", J))
  }
  return(out)
}



# simple regularization of an EM sample 
reg.mc.em.samp <- function(latent, mu){
  n_latent <- length(latent)
  mixture.prop <- mu
  mixture.n <- round(n_latent*mixture.prop) # round to nearest integer 
  uniform.mixture <- runif(mixture.n)
  out <- c(latent,uniform.mixture)
  return(out)
}


# TODO: create a new object which includes the quantiles and the corresponding levels 
likelihood.value <- function(p.hat, tau, latent.trait.quantiles, N, ker, h, J = 100000, mu = 0){
  
  
  # 
  # p.hat = p.hat
  # tau = tau
  # latent.trait.quantiles = tau
  # N = N
  # ker = ker
  # h = h
  # J = 100000
  n_parallel <- 10000
  n_tau <- length(tau) - 1
  weights <- c()
  for(i in 1:n_tau){
    weights[i] <- tau[i + 1] - tau[i]
  }
  
  y.model.samp <- c()
  
  while(length(y.model.samp) < J){
    latent.idx <- sample(1:n_tau, size = n_parallel, replace = TRUE, prob = weights)
    low_bounds <- latent.trait.quantiles[latent.idx]
    high_bounds <- latent.trait.quantiles[latent.idx + 1]
    latent.sample <- runif(n_parallel, min = low_bounds, max = high_bounds)
    
    y.model.samp.tmp <- lapply(latent.sample, function(x){
      i = 0:N
      assumption.weights <- ker((i - N*x)/h)
      out <- sample(0:N, size = 1, replace = TRUE, prob = assumption.weights)
      return(out)
    })
    y.model.samp.tmp <- unlist(y.model.samp.tmp)
    y.model.samp <- c(y.model.samp, y.model.samp.tmp)
  }
  
  like.val <- 0 
  for(i in 0:(N)){
    like.val.tmp <- p.hat[i+1]*log(mean(y.model.samp == i))
    if(is.nan(like.val.tmp)){
      like.val.tmp <- 0
    }
    like.val <- like.val + like.val.tmp
  }
  if(mu > 0){
    
    widths <- c()
    dens.heights <- c()
    for(i in 1:n_tau){
      widths[i] <- (latent.trait.quantiles[i + 1] - latent.trait.quantiles[i])
      dens.heights[i] <- weights[i]/widths[i]
    }
    reg.lik.val <- sum(log(dens.heights)*widths)
    like.val <- like.val + reg.lik.val
  }
  return(like.val)
}


# TODO: create a new object which includes the quantiles and the corresponding levels, this will also be the output of this function
nonpar.em.fit <- function(p.hat, tau, latent.trait.quantiles, N, ker, J, h, mu = 0, verbose = F){
  threshold <- 10^(-6)
  
  
  current.quantiles <- latent.trait.quantiles
  n.quantiles <- length(current.quantiles)
  
  p.ma <- p.ma.from.latent(tau = tau,
                           latent.trait.quantiles = current.quantiles,
                           ker = ker,h = h,N = N,
                           mu = mu,
                           numeric.points = 100)
  
  
  prev.likelihood <- compute.likelihood(p.hat = p.hat, 
                                        p.ma = p.ma, tau = tau, 
                                        latent.trait.quantiles = current.quantiles,
                                        N = N, mu = mu)
  diff.likelihood <- Inf
  em.steps <- 1
  while(diff.likelihood >= threshold){
    latent.sample <- mc.em.samp(p.hat = p.hat, tau = tau, latent.trait.quantiles = current.quantiles, 
                                N = N, ker = ker, J = J, h = h)
    latent.sample <- reg.mc.em.samp(latent.sample, mu)
    current.quantiles <- quantile(latent.sample, tau, names = FALSE)
    current.quantiles[1] <- 0
    current.quantiles[n.quantiles] <- 1
    p.ma <- p.ma.from.latent(tau = tau,
                             latent.trait.quantiles = current.quantiles,
                             ker = ker,h = h,N = N,
                             mu = mu,
                             numeric.points = 100)
    new.likelihood <- compute.likelihood(p.hat = p.hat, 
                                         p.ma = p.ma, tau = tau, 
                                         latent.trait.quantiles = current.quantiles,
                                         N = N, mu = mu.tmp)
    
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




# TODO: unify the computational methods for the A matrix 
A.matrix.compute <- function(R_bins,  N, ker, h, numeric.eval = TRUE, numeric.points = 100, n.mc.samp = 100){
  #R_bins <- 1000
  #n.mc.samp <- 3000 # monte carlo samples for accuracy of A matrix computation 
  A_matrix <- matrix(data = NA, nrow = N + 1, ncol = R_bins)
  if(numeric.eval){
    for(i in 1:R_bins){
      design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
      A.row <- rep(0, N+1)
      for(j in 1:numeric.points){
        y = 0:N
        weights <- ker((y - N*design.points[j])/h)
        A.row <- A.row + weights/(sum(weights))
      }
      A.row <- A.row/sum(A.row)
      A_matrix[,i] = A.row
      cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
    }
    
  }
  else {
    for(i in 1:R_bins){
      uniform.samp <- runif(n.mc.samp, min = (i - 1)/R_bins, max = (i)/R_bins)
      
      y.model.samp <- lapply(uniform.samp, function(x){
        y = 0:N
        assumption.weights <- ker((y - N*x)/h)
        out <- sample(0:N, size = 1, replace = TRUE, prob = assumption.weights)
        return(out)
      })
      y.model.samp <- unlist(y.model.samp)
      for(j in 0:(N)){
        A_matrix[j + 1,i] = mean(y.model.samp == j) 
      }
      cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
    }
  }
  
  return(A_matrix)
}



A.matrix.binomial <- function(R_bins,  N, numeric.points = 100){
  #R_bins <- 1000
  #n.mc.samp <- 3000 # monte carlo samples for accuracy of A matrix computation 
  A_matrix <- matrix(data = NA, nrow = N + 1, ncol = R_bins)
  
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    A.row <- rep(0, N+1)
    for(j in 1:numeric.points){
      y = 0:N
      weights <- dbinom(y, size = N, prob = design.points[j])
      A.row <- A.row + weights/(sum(weights))
    }
    A.row <- A.row/sum(A.row)
    A_matrix[,i] = A.row
    cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
  }
  return(A_matrix)
}



# TODO: document
kl.div <- function(p, q){
  K <- length(p) 
  out <- 0 
  for( k in 1:K){
    out.tmp <- p[k]*log(p[k]/q[k])
    if(!is.nan(out.tmp)){
      out <- out +out.tmp
    }
  }
  out <- out 
  return(out)
}


# TODO: document
tv_norm.dist <- function(p,q){
  (1/2)*sum(abs(p - q))
}


# TODO: Clarify which one of each of these is useful 
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






# TODO: simplify the input and output
# TODO: design separate functions for showing the plot as well as the feasibility test

numeric.latent.fit <- function(p.hat, A.matrix, mu, n.samples, show.plot = FALSE, feasibility.test = FALSE){
  require(CVXR)
  require(ggplot2)
  require(multChernoff)
  
  
  if(missing(n.samples)){
    n.samples <- NA
  }
  R_bins <- ncol(A.matrix)
  N <- nrow(A.matrix) - 1
  theta <- Variable(R_bins, name = "latent discretized distribution")
  obs.dist <- Variable(N + 1, name = "model observed distribution")
  
  data.obj <-  t(p.hat) %*% log(obs.dist)
  #data.obj <- -kl_div(p.hat,obs.dist)
  #pen.obj <- (mu/R_bins)*t(rep(1, R_bins)) %*% log(theta) 
  pen.obj <- (mu/R_bins)*t(rep(1, R_bins)) %*% log(theta) 
  
  constraints <- list(
    obs.dist == A.matrix %*% theta,
    sum(theta) <= 1,
    mu/(R_bins*(1 + mu)) <= theta
  )
  
  obj.arg <- data.obj + pen.obj
  obj <- Maximize(obj.arg)
  p <- Problem(obj, constraints)
  
  
  #value(theta) <- rep(1/R_bins, R_bins) # initial guess of a uniform distribution
  result <- solve(p, solver = "MOSEK")
  #result <- solve(p, verbose = TRUE)
  
  p.m <- result$getValue(theta)
  p.ma <- result$getValue(obs.dist)
  out.list <- list("latent" = p.m, "observed" = p.ma)
  if(show.plot){
    p.max <- max(c(p.hat,p.ma))
    y.max <- p.max*1.15
    plot.ptrue <- ggplot(data = NULL, aes(x = seq(0,N, length.out = N + 1), y = p.hat)) + geom_point() + ggtitle("Empirical Observed Score") + xlab(expression( widehat(p) )) + ylab("Probability Mass") + ylim(0,y.max)
    plot.pma <- ggplot(data = NULL, aes(x = seq(0,N, length.out = N + 1), y = p.ma)) + geom_point() +  ggtitle(paste0("Observed Score Estimate \u00b5 = ", mu)) + xlab(bquote("p"[MA])) + ylab("Probability Mass") + ylim(0,y.max)
    plot.pm <- ggplot(data = NULL, aes(x = seq(0,1, length.out = R_bins), y = R_bins*p.m)) + geom_line() + ggtitle(paste0("Latent Estimate \u00b5 = ", mu)) + xlab("Gamma") + ylab("Density") 
    grid.arrange(plot.pm,                             
                 arrangeGrob(plot.ptrue, plot.pma, ncol = 2), 
                 nrow = 2)   
  } 
  # requires mu = 0 for the feasibility test 
  if(feasibility.test & mu == 0){
    
    lrt.val <- 2*n.samples*kl.div(p.hat, p.ma)
    k <- length(p.hat)
    if(lrt.val == Inf){
      p.feasibility <- 0
    } else {
      p.feasibility <- tailProbBound(x = lrt.val, k = k, n = n.samples)
      if(p.feasibility > 1){
        p.feasibility <- 1
      }
      if(p.feasibility < 0){
        p.feasibility <- 0
      }
    }
    
    out.list <- list("latent" = p.m, "observed" = p.ma, "p_value" = p.feasibility)
  }
  
  
  return(out.list)
  
}




# TODO: rejoin with the other functions 
A.matrix.compute.two.obs <- function(R_bins,  N, ker, h, two.obs.mapping, numeric.points = 100){
  
  A_matrix <- matrix(data = NA, nrow = (N + 1)^2, ncol = R_bins)
  # 3d Array for faster computation. 
  A_3D <- array(NA, dim = c(N+1,N+1,R_bins))
  
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    
    A.block <- array(0, dim = c(N+1,N+1))
    y1 = 0:N
    y2 = 0:N
    
    for(j in 1:numeric.points){
      
      weights <- outer(ker((y1 - N*design.points[j])/h), ker((y2 - N*design.points[j])/h), "*")
      A.block <- A.block + weights/(sum(weights))
    }
    A_3D[,,i] <- A.block/(sum(A.block))
    cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
  }
  
  for(k in 1:((N+1)^2)){
    
    y1.idx = two.obs.mapping[k,1]
    y2.idx = two.obs.mapping[k,2]
    
    
    A_matrix[k,] <- A_3D[y1.idx + 1,y2.idx + 1,]
  }
  
  return(A_matrix)
}


# TODO: unify with the general case
A.matrix.binomial.two.obs <- function(R_bins, N,two.obs.mapping, numeric.points = 100){
  
  A_matrix <- matrix(data = NA, nrow = (N + 1)^2, ncol = R_bins)
  # 3d Array for faster computation. 
  A_3D <- array(NA, dim = c(N+1,N+1,R_bins))
  
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    
    A.block <- array(0, dim = c(N+1,N+1))
    y1 = 0:N
    y2 = 0:N
    
    for(j in 1:numeric.points){
      
      weights <- outer(dbinom(y1, size = N, prob = design.points[j]), dbinom(y2, size = N, prob = design.points[j]), "*")
      
      A.block <- A.block + weights/(sum(weights))
    }
   
    A_3D[,,i] <- A.block/(sum(A.block))
    cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
  }
  
  for(k in 1:((N+1)^2)){
    
    y1.idx = two.obs.mapping[k,1]
    y2.idx = two.obs.mapping[k,2]
    
    
    A_matrix[k,] <- A_3D[y1.idx + 1,y2.idx + 1,]
  }
  
  return(A_matrix)
}


# TODO: unify with the cases, leave tensor or matrix as options for the outputs
A.two.obs.tensor.compute <- function(R_bins,  N, ker, h, numeric.points = 100){
  
  # 3d Array for faster computation. 
  A_3D <- array(NA, dim = c(N+1,N+1,R_bins))
  
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    
    A.block <- array(0, dim = c(N+1,N+1))
    y1 = 0:N
    y2 = 0:N
    
    for(j in 1:numeric.points){
      
      weights <- outer(ker((y1 - N*design.points[j])/h), ker((y2 - N*design.points[j])/h), "*")
      A.block <- A.block + weights/(sum(weights))
    }
    A_3D[,,i] <- A.block/(sum(A.block))
    cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
  }
  
  return(A_3D)
}


A.two.obs.tensor.binomial <- function(R_bins,  N, numeric.points = 100){
  
  # 3d Array for faster computation. 
  A_3D <- array(NA, dim = c(N+1,N+1,R_bins))
  
  for(i in 1:R_bins){
    design.points <- seq((i - 1)/R_bins, (i)/R_bins, length.out = numeric.points)
    
    A.block <- array(0, dim = c(N+1,N+1))
    y1 = 0:N
    y2 = 0:N
    
    for(j in 1:numeric.points){
      
      weights <- outer(dbinom(y1, size = N, prob = design.points[j]), dbinom(y2, size = N, prob = design.points[j]), "*")
      A.block <- A.block + weights/(sum(weights))
    }
    A_3D[,,i] <- A.block/(sum(A.block))
    cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
  }
  
  return(A_3D)
}


two.samp.log.likelihood <- function(y.pair, A.two.sample.tensor, latent.mixture.list){
  M <- nrow(y.pair)
  log.likelihood <- 0
  for(m in 1:M){
    like.tmp <- t(A.two.sample.tensor[y.pair[m,1] +1,y.pair[m,2] + 1,]) %*% latent.mixture.list[[m]]
    log.likelihood <- log.likelihood + log(like.tmp)
  }
  return(log.likelihood)
}

# TO DO: simplify this method.  I.e. create a categorical variable which covers all the age ranges.  
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


# TODO: unify with the general case
tau.set <- function(latent.mixture){
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


# TODO: rename, correct characters
# TODO: unify with the binomial
intrinsic.variability.samp <- function(pair.obs, latent.mixture, n.samp, N, ker,h, p.hat){

  
  tau <- tau.set(latent.mixture)
  latent.quantiles <- seq(0,1, length.out = length(tau))
  
  if(missing(p.hat)){
    p.hat <- rep(1,N + 1)
    p.hat <- p.hat/sum(p.hat)
  }
  
  latent.samp <- conditional.samp(y.cond = pair.obs[1], n.samp = n.samp,
                                  p.hat = p.hat, tau = tau, latent.trait.quantiles = latent.quantiles,
                                  N = N, ker = ker, h = h)
  
  
  y.model.samp <- lapply(latent.samp, function(x){
    i = 0:N
    assumption.weights <- ker((i - N*x)/h)
    out <- sample(0:N, size = 1, replace = TRUE, prob = assumption.weights)
    return(out)
  })
  y.model.samp <- unlist(y.model.samp)
  
  d.true <- pair.obs[2] - pair.obs[1]
  d.sim <- y.model.samp - pair.obs[1]
  
  out.list <- list("true.diff" = d.true, "sim.diff" = d.sim)
  return(out.list)
}



# TODO: simplify to a single version 

intrinsic.variability <- function(y.true.frame, latent.mix.list, model.observed.list, n.samp, N, ker, h, show.plot = FALSE, parallel = TRUE){
  d.true <- y.true.frame[,2] - y.true.frame[,1]
  
  d.sim <- c()
  
  
  if(parallel){
    idx <- 1:nrow(y.true.frame)
    d.sim.list <- lapply(idx, function(x){
      latent.mix <- latent.mix.list[[x]]
      train.p.hat <- model.observed.list[[x]]
      intrinsic.samp <- intrinsic.variability.samp(pair.obs = as.numeric(y.true.frame[x,]), latent.mixture = latent.mix,
                                                   n.samp = n.samp, N = N, ker = ker, h = h, p.hat = train.p.hat)
      
      return(intrinsic.samp$sim.diff)
    })
    d.sim <- unlist(d.sim.list)
  } else {
    for(i in 1:nrow(y.true.frame)){
      latent.mix <- latent.mix.list[[i]]
      train.p.hat <- model.observed.list[[i]]
      intrinsic.samp <- intrinsic.variability.samp(pair.obs = as.numeric(y.true.frame[i,]), latent.mixture = latent.mix,
                                                   n.samp = n.samp, N = N, ker = ker, h = h, p.hat = train.p.hat)
      
      d.true <- c(d.true,intrinsic.samp$true.diff)
      d.sim <- c(d.sim,intrinsic.samp$sim.diff)
      cat(paste0("Latent Sample: ", i,"/",nrow(y.true.frame)), end="\r")
    }
  }
  if(show.plot){
    iv.plot.true <- ggplot(data = NULL, aes(x = d.true)) + geom_histogram(bins = 2*N + 1) + ggtitle("Empirical Variability") + xlab("Score Difference") + ylab("Probability Mass")
    iv.plot.sim <- ggplot(data = NULL, aes(x = d.sim)) + geom_histogram(bins = 2*N + 1) + ggtitle("Model Variability") + xlab("Score Difference") + ylab("Probability Mass")
    grid.arrange(arrangeGrob(iv.plot.true, iv.plot.sim, ncol = 2),
                 nrow = 1)
  }
  
  
  
  int.norm <- tv_norm(d.true,d.sim)
  return(int.norm)
}






# TODO: should not require a different function for this task 
intrinsic.variability.binom.samp <- function(pair.obs, latent.mixture, n.samp, N, p.hat){
  
  
  tau <- tau.set(latent.mixture)
  latent.quantiles <- seq(0,1, length.out = length(tau))
  
  if(missing(p.hat)){
    p.hat <- rep(1,N + 1)
    p.hat <- p.hat/sum(p.hat)
  }
  
  latent.samp <- conditional.binom.samp(y.cond = pair.obs[1], n.samp = n.samp,
                                  p.hat = p.hat, tau = tau, latent.trait.quantiles = latent.quantiles,
                                  N = N)
  
  
  y.model.samp <- lapply(latent.samp, function(x){
    out <- rbinom(n = 1, size = N, prob = x)
    return(out)
  })
  y.model.samp <- unlist(y.model.samp)
  
  d.true <- pair.obs[2] - pair.obs[1]
  d.sim <- y.model.samp - pair.obs[1]
  
  out.list <- list("true.diff" = d.true, "sim.diff" = d.sim)
  return(out.list)
}



# TODO: should not require a different function for this task 
intrinsic.variability.binom <- function(y.true.frame, latent.mix.list, model.observed.list, n.samp, N,show.plot = FALSE, parallel = TRUE){
  d.true <- y.true.frame[,2] - y.true.frame[,1]
  
  d.sim <- c()
  
  
  if(parallel){
    idx <- 1:nrow(y.true.frame)
    d.sim.list <- lapply(idx, function(x){
      latent.mix <- latent.mix.list[[x]]
      train.p.hat <- model.observed.list[[x]]
      intrinsic.samp <- intrinsic.variability.binom.samp(pair.obs = as.numeric(y.true.frame[x,]), latent.mixture = latent.mix,
                                                   n.samp = n.samp, N = N, p.hat = train.p.hat)
      
      return(intrinsic.samp$sim.diff)
    })
    d.sim <- unlist(d.sim.list)
  } else {
    for(i in 1:nrow(y.true.frame)){
      latent.mix <- latent.mix.list[[i]]
      train.p.hat <- model.observed.list[[i]]
      intrinsic.samp <- intrinsic.variability.binom.samp(pair.obs = as.numeric(y.true.frame[i,]), latent.mixture = latent.mix,
                                                         n.samp = n.samp, N = N, p.hat = train.p.hat)
      
      d.true <- c(d.true,intrinsic.samp$true.diff)
      d.sim <- c(d.sim,intrinsic.samp$sim.diff)
      cat(paste0("Latent Sample: ", i,"/",nrow(y.true.frame)), end="\r")
    }
  }
  if(show.plot){
    iv.plot.true <- ggplot(data = NULL, aes(x = d.true)) + geom_histogram(bins = 2*N + 1) + ggtitle("Empirical Variability") + xlab("Score Difference") + ylab("Probability Mass")
    iv.plot.sim <- ggplot(data = NULL, aes(x = d.sim)) + geom_histogram(bins = 2*N + 1) + ggtitle("Model Variability") + xlab("Score Difference") + ylab("Probability Mass")
    grid.arrange(arrangeGrob(iv.plot.true, iv.plot.sim, ncol = 2),
                 nrow = 1)
  }
  
  
  
  int.norm <- tv_norm(d.true,d.sim)
  return(int.norm)
}


# TODO: document
placeholder.est <-function(){
  out <- list('latent' = NULL, 'observed' = NULL, 'p_value' = 0)
  return(out)
}


# TODO: document
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

# TODO: document
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

## Outputs the weight in a discretized bin region 
latent.mix.val <- function(x,latent.mixture){
  R_bins <- length(latent.mixture)
  binpoints <- seq(0,1, length.out = R_bins + 1)
  
  dens.height <- sapply(x, function(z){
    idx <- which.min(binpoints <= z) - 1 # Selects proper bin for density height 
    out <- latent.mixture[idx]*R_bins
    if(z == 1){
      out <- latent.mixture[R_bins]
    } else if (z < 0 | z > 1){
      out <- NA
    }
    return(out)})
  
  return(dens.height)
}


# TODO: verify and delete 
wasserstein1 <- function(z.true, z.pred){
  abs.diff <- abs(z.true - z.pred)
  out <- mean(abs.diff)
  return(out)
}


# TODO: document and come up with a consistent naming scheme
# TODO: replace . with _ 
# TODO: document, converting to numeric conditional distribution
conversion.numeric <- function(y, z, latent.mixture.y, latent.mixture.z, 
                               Ny, Nz, ker.y, ker.z, hy, hz, p.ma.y, 
                               binomial.y = F, binomial.z = F){
  
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
  
  tau.y <- tau.set(latent.mixture.y)
  latent.quantiles.y <- seq(0,1, length.out = length(tau.y))
  
  tau.z <- tau.set(latent.mixture.z)
  latent.quantiles.z <- seq(0,1, length.out = length(tau.z))
  
  R_bins <- length(latent.mixture.y)
  latent.points <- sapply(1:R_bins, function(z){
    gam <- (latent.quantiles.y[z] + latent.quantiles.y[z+1])/2
    return(gam)
  })
  
  ######################################################
  ######################################################
  
  if(!binomial.y){
    p.ay <- function(y,gam){
      i <- 0:Ny
      out <- ker.y((y - Ny*gam)/hy)/sum(ker.y((i - Ny*gam)/hy))
    }
  }
  
  
  if(!binomial.z){
    p.az <- function(z,zeta){
      j <- 0:Nz
      out <- ker.z((z - Nz*zeta)/hz)/sum(ker.z((j - Nz*zeta)/hz))
    }
  }
  
  if(binomial.y){
    p.ay <- function(y,gam){
      out <- dbinom(x = y, size = Ny, prob = gam)
    }
  }
  
  
  if(binomial.z){
    p.az <- function(z,zeta){
      out <- dbinom(x = z, size = Nz, prob = zeta)
    }
  }
  
  
  joint.dist.grid  <- sapply(1:R_bins, function(gam.idx){
    gam <- latent.points[gam.idx]
    weight <- latent.mixture.y[gam.idx]
    zeta.point <- qf_L(cdf_L(gam, latent.quantiles.y, tau.y), latent.quantiles.z, tau.z)
    
    out.y <- p.ay(y,gam)
    out.z <- p.az(z,zeta.point)
    out <- out.y*out.z*weight
    return(out)
  })
  denom <- p.ma.y[y + 1]
  cond.out <- sum(joint.dist.grid)/(denom)
  
  return(cond.out)
}


# TODO: verify if this is ever even necessary 

conversion.samp <- function(y, z, latent.mixture.y, latent.mixture.z, n.samp, Ny, Nz, ker.y, ker.z, hy, hz, p.hat.y){
  tau.y <- tau.set(latent.mixture.y)
  latent.quantiles.y <- seq(0,1, length.out = length(tau.y))
  
  tau.z <- tau.set(latent.mixture.z)
  latent.quantiles.z <- seq(0,1, length.out = length(tau.z))
  
  
  if(missing(p.hat.y)){
    p.hat <- rep(1,Ny + 1)
    p.hat <- p.hat/sum(p.hat)
  }
  
  # sample latent gamma (y's latent variable)
  latent.samp.y <- conditional.samp(y.cond = y, n.samp = n.samp,
                                    p.hat = p.hat.y, 
                                    tau = tau.y, 
                                    latent.trait.quantiles = latent.quantiles.y,
                                    N = Ny, ker = ker.y, h = hy)
  
  
  # converting scores by the linear interpolation of quantiles
  
  latent.samp.z <- qf_L(cdf_L(latent.samp.y, latent.quantiles.y, tau.y), latent.quantiles.z, tau.z)
  
  z.model.samp <- sapply(latent.samp.z, function(x){
    i = 0:Nz
    assumption.weights <- ker.z((i - Nz*x)/hz)
    out <- sample(0:Nz, size = 1, replace = TRUE, prob = assumption.weights)
    return(out)
  })
  
  
  diff.true <- z - y
  diff.sim <- z.model.samp - y
  
  z.true <- z
  z.sim <- z.model.samp
  
  out.list <- list("true.diff" = diff.true, "sim.diff" = diff.sim, "z.true" = z, "z.sim" = z.model.samp)
  return(out.list)
}


# TODO: simplify greatly, we only care about the cross entropy
# TODO: rename


score.conversion <- function(test.pairs, latent.mix.list.y, model.observed.list.y,
                             latent.mix.list.z, model.observed.list.z, Ny, ker.y, hy,
                             Nz, ker.z, hz, joint.prob,n.samp = 5, method = "CrossEntropy", 
                             binomial.y = F, binomial.z = F){
  # output error if binomial false and missing kernel and bandwidth
  
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
  
  if(missing(joint.prob)){
    joint.prob <- matrix(data = 1, nrow = Ny + 1, ncol = Nz + 1)
  }
  d.true <- test.pairs[,2] - test.pairs[,1]
  z.true <- test.pairs[,2]
  d.sim <- c()
  z.sim <- c()
  
  
  idx <- 1:nrow(test.pairs)
  res <- sapply(idx, function(x){
    latent.mix.y <- latent.mix.list.y[[x]]
    p.ma.y <- model.observed.list.y[[x]]
    
    latent.mix.z <- latent.mix.list.z[[x]]
    p.ma.z <- model.observed.list.z[[x]]
    
    if(method == "CrossEntropy"){
      
      if(!binomial.y & !binomial.z){
        conv.prob <- conversion.numeric(y = test.pairs[x,1], z = test.pairs[x,2],
                                        latent.mixture.y = latent.mix.y,
                                        latent.mixture.z = latent.mix.z,
                                        Ny = Ny, Nz = Nz,
                                        ker.y = ker.y, ker.z = ker.z,
                                        hy = hy, hz = hz, p.ma.y = p.ma.y)
      }
      
      if(!binomial.y & binomial.z){
        conv.prob <- conversion.numeric(y = test.pairs[x,1], z = test.pairs[x,2],
                                        latent.mixture.y = latent.mix.y,
                                        latent.mixture.z = latent.mix.z,
                                        Ny = Ny, Nz = Nz,
                                        ker.y = ker.y, 
                                        hy = hy, p.ma.y = p.ma.y, binomial.z = T)
      }
      
      if(binomial.y & !binomial.z){
        conv.prob <- conversion.numeric(y = test.pairs[x,1], z = test.pairs[x,2],
                                        latent.mixture.y = latent.mix.y,
                                        latent.mixture.z = latent.mix.z,
                                        Ny = Ny, Nz = Nz,
                                        ker.z = ker.z,
                                        hz = hz, p.ma.y = p.ma.y, binomial.y = T)
      }
      
      if(binomial.y & binomial.z){
        conv.prob <- conversion.numeric(y = test.pairs[x,1], z = test.pairs[x,2],
                                        latent.mixture.y = latent.mix.y,
                                        latent.mixture.z = latent.mix.z,
                                        Ny = Ny, Nz = Nz,
                                        p.ma.y = p.ma.y, binomial.y = T, binomial.z = T)
      }
      
      joint.weight <- joint.prob[test.pairs[x,1] + 1, test.pairs[x,2] + 1]
      out <- -joint.weight*log(conv.prob)
      
    } else if (method == "ExtrinsicVariability"){
      conversions <- conversion.samp(y = test.pairs[x,1], z = test.pairs[x,2],
                                     latent.mixture.y = latent.mix.y,
                                     latent.mixture.z = latent.mix.z,
                                     n.samp = n.samp, Ny = Ny, Nz = Nz,
                                     ker.y = ker.y, ker.z = ker.z,
                                     hy = hy, hz = hz, p.hat.y = p.ma.y)
      out <- conversions$sim.diff
    } else if (method == "Wasserstein1"){
      conversions <- conversion.samp(y = test.pairs[x,1], z = test.pairs[x,2],
                                     latent.mixture.y = latent.mix.y,
                                     latent.mixture.z = latent.mix.z,
                                     n.samp = n.samp, Ny = Ny, Nz = Nz,
                                     ker.y = ker.y, ker.z = ker.z,
                                     hy = hy, hz = hz, p.hat.y = p.ma.y)
      pred.samp <- conversions$z.sim
      z.true <- conversions$z.true
      out <- wasserstein1(z.true, pred.samp)
    }
    cat(paste0("Test sample ", x, '/',nrow(test.pairs)), end = '\r')
    return(out)
  })
  
  if(method == "CrossEntropy" ){
    out <- sum(res)
  } else if (method == "ExtrinsicVariability"){
    out <- tv_norm(d.true,res)
  } else if (method == "Wasserstein1"){
    out <- sum(res)
  }
  
  cat(end = '\n')
  return(out)
}



# TODO: rename,  

p.ma.from.latent <- function(tau,latent.trait.quantiles,
                             ker,h,N,
                             mu = 0,
                             numeric.points = 100){
  
  
  
  
  
  R_bins <- length(latent.trait.quantiles) - 1
  A_matrix <- matrix(data = NA, nrow = N + 1, ncol = R_bins)
  for(i in 1:R_bins){
    design.points <- seq(latent.trait.quantiles[i], 
                         latent.trait.quantiles[i+1], 
                         length.out = numeric.points)
    A.row <- rep(0, N+1)
    for(j in 1:numeric.points){
      y = 0:N
      weights <- ker((y - N*design.points[j])/h)
      A.row <- A.row + weights/(sum(weights))
    }
    A.row <- A.row/sum(A.row)
    A_matrix[,i] = A.row
    #cat(paste0("A Matrix Computed Row: ", i,"/",R_bins), end="\r")
  }
  
  weights <- c()
  for(i in 1:R_bins){
    weights[i] <- tau[i + 1] - tau[i]
  }
  p.ma <- A_matrix %*% weights
  return(p.ma)
}


# TODO: verify if this is ever even necessary 

quantile.conversion.prob <- function(y,quantile_map_yz, decay_rate, Nz, n.samp = 100000){
  z.pred <- quantile_map_yz[y + 1]
  
  z.scale <- 0:Nz
  probs <- exp(-decay_rate*abs(z.pred - z.scale))
  
  out <- probs/sum(probs)
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



# TODO: rename, to psi function 
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
