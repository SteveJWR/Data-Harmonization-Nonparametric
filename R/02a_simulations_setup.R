
source("R/01_functions.R")

# picking a set of grid tuning parameters for simulations. 
h.set<-  c(0.25,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0)  #c(1,2) #c(0.25,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 

ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)# list(gaussian_kernel, epanechnikov_kernel) #list(gaussian_kernel, exponential_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 3))) #c(0,0.1) #c(0,exp(seq(log(0.001),log(0.1), length.out = 3)))
mu.set.long <- c(0,exp(seq(log(0.001),log(.1), length.out = 15))) # For selection of Mu 

# Find these to be optimal
mu.y <- 0.01930
mu.z <- 0.01390 

mu.set.conversion <-  c(0,0.1)# c(0,exp(seq(log(0.001),log(mu.z), length.out = 3)), exp(seq(log(mu.y),log(.3), length.out = 3)))

h.set.feasibility.test <- c(0.25,0.5,1,2,4,7,10,14,20,25,30)

h.set.speed.test <- c(0.25,0.5, 0.75,1.0,2.0 ,3.0, 5.0, 8.0) 

# number of individuals in the data set
dataset.size <- 100
# number of simulations 
n.sims <- 100
n.sims.feasibility.test <- 300

# list of possible tuning parameters for the intrinsic variability 
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))

# compact support kernels must have h >= 1
idx1 <- hyper.param.idx[,1]  %in% which(h.set <= 1)
idx2 <- hyper.param.idx[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.idx <- hyper.param.idx[!(idx1 & idx2),]


# list of possible tuning parameters for the conversion
hyper.param.conversion.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set.conversion))


# compact support kernels must have h > 1
idx1 <- hyper.param.conversion.idx[,1]  %in% which(h.set <= 1)
idx2 <- hyper.param.conversion.idx[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.conversion.idx <- hyper.param.conversion.idx[!(idx1 & idx2),]

# Indexing guide for hyper parameter with NPMLE i.e. mu = 0
hyper.param.conversion.ml.idx <- hyper.param.conversion.idx[hyper.param.conversion.idx[,3] == 1,]


N = 30 
Ny = N 
Nz = N
# set of possible matrices corresponding to kernel bandwidth pairs
A.matrix.set <- readRDS("Data/A_matrix_set.RDS")
if(!exists("A.matrix.set")){
  # Pre computing A matrices
  A.matrix.set <- list()
  for(j in 1:length(h.set)){
    h.tmp <- h.set[j]
    A.matrix.set.tmp <- list()
    for(k in 1:length(ker.set)){
      
      ker.tmp <- ker.set[[k]]
      cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
      A.mat <- compute_A_matrix(R_bins = 1000, cond = cond.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.mat[is.nan(A.mat)] <- 0
      A.matrix.set.tmp[[k]] <- A.mat 
      
    }
    A.matrix.set[[j]] <- A.matrix.set.tmp
  }
}
#saveRDS(A.matrix.set, "Data/A_matrix_set.RDS")


# Simulation 1 Intrinsic Variability setup 
# ---------------------------


### true model parameters (model 1)
k.model1 <- gaussian_kernel
h.model1 <- 2

beta1.model1 <- 12 
beta2.model1 <- 5

cond.model1 <- conditional_mkm(N, ker = k.model1, h = h.model1)
A.matrix.model1 <- compute_A_matrix(R_bins = 1000, cond = cond.model1, numeric.points = 400)
A.tensor.model1 <- compute_A_two_obs_tensor(R_bins = 1000, cond = cond.model1, numeric.points = 400)







# Simulation 2 Intrinsic Variability setup 
# ---------------------------


### true model parameters (model 1)
k.model2 <- exponential_kernel
h.model2 <- 1

beta1.model2 <- 6 
beta2.model2 <- 6

cond.model2 <- conditional_mkm(N, ker = k.model2, h = h.model2)
A.matrix.model2 <- compute_A_matrix(R_bins = 1000, cond = cond.model2, numeric.points = 400)
A.tensor.model2 <- compute_A_two_obs_tensor(R_bins = 1000, cond = cond.model2, numeric.points = 400)





# ---- Feasibility Test Grid options ---- 

# set of possible matrices corresponding to kernel bandwidth pairs
A.matrix.set.feasibility <- readRDS("Data/A_matrix_set_feasibility.RDS")
if(!exists("A.matrix.set.feasibility")){
  # Pre computing A matrices
  A.matrix.set.feasibility <- list()
  for(j in 1:length(h.set)){
    h.tmp <- h.set.feasibility[j]
    A.matrix.set.feasibility.tmp <- list()
    for(k in 1:length(ker.set)){
      
      ker.tmp <- ker.set[[k]]
      cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
      A.mat <- compute_A_matrix(R_bins = 1000, cond = cond.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.mat[is.nan(A.mat)] <- 0
      A.matrix.set.tmp[[k]] <- A.mat 
      
    }
    A.matrix.set[[j]] <- A.matrix.set.tmp
  }
}
#saveRDS(A.matrix.set.feasibility, "Data/A_matrix_set_feasibility.RDS")




# set of possible tensors corresponding to kernel bandwidth pairs
A.tensor.set.feasibility <- readRDS("Data/A_tensor_set_feasibility.RDS")
if(!exists("A.tensor.set.feasibility")){
  # Pre computing A tensors
  A.tensor.set.feasibility <- list()
  for(j in 1:length(h.set)){
    h.tmp <- h.set.feasibility[j]
    A.tensor.set.feasibility.tmp <- list()
    for(k in 1:length(ker.set)){
      
      ker.tmp <- ker.set[[k]]
      cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
      A.tensor <- compute_A_two_obs_tensor(R_bins = 1000, cond = cond.tmp, numeric.points = 400)
      
      # replacing NAN with 0 due to numerical rounding error 
      A.tensor[is.nan(A.tensor)] <- 0
      A.tensor.set.tmp[[k]] <- A.tensor 
      
    }
    A.tensor.set[[j]] <- A.tensor.set.tmp
  }
}
#saveRDS(A.tensor.set.feasibility, "Data/A_tensor_set_feasibility.RDS")








## ---- Simulation functions ---- 


simulation_intrinsic_variability_model1 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  intrinsic.variability.each.model <- sapply(1:nrow(hyper.param.idx), function(x){
    
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    l = hyper.param.idx[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- mu.set[l]
    A.matrix <- A.matrix.set[[j]][[k]]
    
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix, mu = mu.tmp)
    
    cond <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    model.observed.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
      model.observed.list[[m]] <- model.estimate$observed
    }
    # intrinsic variability results
    iv.res <- intrinsic_variability(pair.obs = sim.data, latent.mix.list = latent.mix.list, 
                                    model.observed.list = model.observed.list, n.samp = 25, cond = cond)
    
    return(iv.res)
  })
  return(intrinsic.variability.each.model)
}




simulation_intrinsic_variability_model2 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model2, beta2.y = beta2.model2, 
                            cond.y = cond.model2, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  intrinsic.variability.each.model <- sapply(1:nrow(hyper.param.idx), function(x){
    
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    l = hyper.param.idx[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- mu.set[l]
    A.matrix <- A.matrix.set[[j]][[k]]
    
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix, mu = mu.tmp)
    
    cond <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    model.observed.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
      model.observed.list[[m]] <- model.estimate$observed
    }
    # intrinsic variability results
    iv.res <- intrinsic_variability(pair.obs = sim.data, latent.mix.list = latent.mix.list, 
                                    model.observed.list = model.observed.list, n.samp = 25, cond = cond)
    
    return(iv.res)
  })
  return(intrinsic.variability.each.model)
}





## -------- mu selection --------
## --------    Model 1   --------


simulation_smoothing_selection_model1 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  neg.log.like.each.mu <- sapply(mu.set.long, function(x){
    mu.tmp <- x
    
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix.model1, mu = mu.tmp)
    
    cond <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
    }
    # intrinsic variability results
    log.like <- compute_two_obs_loglikelihood(sim.data, A.tensor.model1, latent.mix.list)
    
    
    neg.log.like <- -log.like
    return(neg.log.like)
  })
  return(neg.log.like.each.mu)
}


## --------    Model 2   --------

simulation_smoothing_selection_model2 <- function(sim.number){
  set.seed(sim.number)
  
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model2, beta2.y = beta2.model2, 
                            cond.y = cond.model2, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  
  ## Computes the intrinsic variability across all possible models in our set 
  neg.log.like.each.mu <- sapply(mu.set.long, function(x){
    mu.tmp <- x
    
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix.model2, mu = mu.tmp)
    
    cond <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
    }
    # intrinsic variability results
    log.like <- compute_two_obs_loglikelihood(sim.data, A.tensor.model2, latent.mix.list)
    
    
    neg.log.like <- -log.like
    return(neg.log.like)
  })
  return(neg.log.like.each.mu)
}



## --------   Conversion Cross Entropy   --------

### population test set 


test.grid <- expand.grid(0:Ny, 0:Nz)
p.yz <- compute_joint_dist_beta(beta1.model.y = beta1.model1, beta2.model.y = beta2.model1,
                                beta1.model.z = beta1.model2, beta2.model.z = beta2.model2,
                                cond.y = cond.model1, cond.z = cond.model2, grid.size = 10000)




simulation_conversion_cross_entropy <- function(sim.number){
  set.seed(sim.number)
  # setup of parameters of the simulation
  
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1,
                            beta1.z = beta1.model2, beta2.z = beta2.model2, 
                            cond.z = cond.model2,
                            pair.obs = T)
  
  # only keep the single observations 
  sim.data <- sim.data[,c(1,3)]
  train.p.hat.y <- compute_edf(sim.data[,1], N)
  train.p.hat.z <- compute_edf(sim.data[,2], N)
  
  
  ### Fixed Model Estimate of y under correctly learned model 
  model.estimate.y <- estimate_mixing_numeric(p.hat = train.p.hat.y, A.matrix = A.matrix.model1, mu = mu.y)
  model.estimate.y.ml <- estimate_mixing_numeric(p.hat = train.p.hat.y, A.matrix = A.matrix.model1, mu = 0)
  
  
  ## Computes the intrinsic variability across all possible models in our set 
  conversion.each.model <- sapply(1:nrow(hyper.param.conversion.idx), function(x){
    
    j = hyper.param.conversion.idx[x,1]
    k = hyper.param.conversion.idx[x,2]
    l = hyper.param.conversion.idx[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- mu.set[l]
    A.matrix.z.tmp <- A.matrix.set[[j]][[k]]
    
    model.estimate.z <- estimate_mixing_numeric(p.hat = train.p.hat.z, 
                                                A.matrix = A.matrix.z.tmp, 
                                                mu = mu.tmp)
    
    cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    
    latent.mix.list.y <- list()
    model.observed.list.y <- list()
    latent.mix.list.z <- list()
    model.observed.list.z <- list()
    for (m in 1:nrow(test.grid)) {
      latent.mix.list.y[[m]] <- model.estimate.y$latent
      latent.mix.list.z[[m]] <- model.estimate.z$latent
    }
    
    # conversion population cross entropy 
    ce.pop <- convert_score_metric(test.grid, 
                                   latent.mix.list.y, 
                                   latent.mix.list.z, 
                                   cond.y, 
                                   cond.z, 
                                   joint.prob = p.yz,
                                   grid.size = 200) # TODO: Replace with larger value
    
    return(ce.pop)
  })
  
  
  conversion.each.model.ml <- sapply(1:nrow(hyper.param.conversion.ml.idx), function(x){
    
    j = hyper.param.conversion.ml.idx[x,1]
    k = hyper.param.conversion.ml.idx[x,2]
    l = hyper.param.conversion.ml.idx[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    
    A.matrix.z.tmp <- A.matrix.set[[j]][[k]]
    
    model.estimate.z.ml <- estimate_mixing_numeric(p.hat = train.p.hat.z, 
                                                A.matrix = A.matrix.z.tmp, 
                                                mu = 0)
    
    cond.tmp <- conditional_mkm(N,ker = ker.tmp, h = h.tmp)
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    
    latent.mix.list.y.ml <- list()
    latent.mix.list.z.ml <- list()
    for (m in 1:nrow(test.grid)) {
      latent.mix.list.y.ml[[m]] <- model.estimate.y.ml$latent
      latent.mix.list.z.ml[[m]] <- model.estimate.z.ml$latent
    }
    
    # conversion population cross entropy 
    ce.pop.ml <- convert_score_metric(test.grid, 
                                      latent.mix.list.y.ml, 
                                      latent.mix.list.z.ml, 
                                      cond.y, 
                                      cond.z, 
                                      joint.prob = p.yz,
                                      grid.size = 200) # TODO: Replace with larger value
    
    return(ce.pop.ml)
  })
  
  results <- list("smoothed" = conversion.each.model,
                  "npmle" = conversion.each.model.ml)
  return(results)
}





## feasibility test simulations


simulation_feasibility_tests <- function(sim.number){
  set.seed(sim.number)
  
  
  
  
}











