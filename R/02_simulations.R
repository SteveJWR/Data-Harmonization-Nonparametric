



source("R/DataHarmonizationFunctions.R")

# picking a set of grid tuning parameters for simulations. 
h.set<- c(0.2,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 
h.set.feasibility <- c(0.2,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 
ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 3)))
mu.set.long <- c(0,exp(seq(log(0.001),log(.1), length.out = 15))) # For selection of Mu 

# number of individuals in the data set
dataset.size <- 100
# number of simulations 
n.sims <- 100

N = 30 

#
# Pre computing A matrices
A.matrix.set <- list()
for(j in 1:length(h.set)){
  h.tmp <- h.set[j]
  A.matrix.set.tmp <- list()
  for(k in 1:length(ker.set)){
    
    ker.tmp <- ker.set[[k]]
    
    A.mat <- A.matrix.compute(R_bins = 1000, N = 30, ker = ker.tmp, h = h.tmp, numeric.points = 400)
    A.matrix.set.tmp[[k]] <- A.mat 
    
  }
  A.matrix.set[[j]] <- A.matrix.set.tmp
}


# Simulation 1 Intrinsic Variability 
# ---------------------------


### true model parameters (model 1)
k.model1 <- gaussian_kernel
h.model1 <- 2

beta1.model1 <- 12 
beta2.model1 <- 5

cond.model1 <- conditional_mkm(N, ker = k.model1, h = h.model1)


results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))
# list of possible tuning parameters 
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))

# compact support kernels must have h >= 1
idx1 <- hyper.param.idx[,1]  %in% which(h.set < 1)
idx2 <- hyper.param.idx[,2]  %in% c(3,4)

# Indexing guide for hyper parameter options compared to the sets used 
hyper.param.idx <- hyper.param.idx[!(idx1 & idx2),]


results.list <- sapply(1:n.sims, function(i){
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  train.p.hat <- compute_edf(sim.data[,1], N)
  print(paste0("Dataset Sample: ", i, " of ", n.sims))
  
  
  ## Computes the intrinsic variability across all possible models
  intrinsic.variability.each.model <- sapply(1:nrow(hyper.param.idx), function(x){
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    l = hyper.param.idx[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- mu.set[l]
    A.matrix <- A.matrix.set[[j]][[k]]
    model.estimate <- estimate_mixing_numeric(p.hat = train.p.hat, A.matrix = A.matrix, mu = mu.tmp)
    
    
    # lists of latent and model implied observed distributions for allowing 
    # different latent distributions for each point
    latent.mix.list <- list()
    model.observed.list <- list()
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
      model.observed.list[[m]] <- model.estimate$observed
    }
    
    res <- intrinsic.variability(y.true.frame = sim.data, latent.mix.list = latent.mix.list, 
                                 model.observed.list = model.observed.list, n.samp = 5, N = N,
                                 ker = ker.tmp, h = h.tmp, parallel = T, show.plot = F)
    
    return(res)
  })
  res.set <- unlist(res.list)
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    l = hyper.param.idx[q,3]
    results.array[i,j,k,l] <- res.set[q]
  }
  
  
  
  
})

for(i in 1:n.sims){
  sim.data <- simulate_beta(n.ind = dataset.size, n.obs.per.ind = 2, 
                            beta1.y = beta1.model1, beta2.y = beta2.model1, 
                            cond.y = cond.model1, pair.obs = F)
  
  
  
  train.p.hat <- compute_edf(sim.data[,1], N)
  print(paste0("Dataset Sample: ", i, " of ", n.sims))
  
  ## Computes the intrinsic variability across all possible models
  intrins.var.all.models <- sapply(1:nrow(hyper.param.idx), function(x){
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    l = hyper.param.idx[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- mu.set[l]
    A.matrix <- A.matrix.set[[j]][[k]]
    model.estimate <- numeric.latent.fit(p.hat = train.p.hat, A.matrix = A.matrix, mu = mu.tmp, show.plot = F)
    
    latent.mix.list <- list()
    model.observed.list <- list()
    
    for (m in 1:dataset.size) {
      latent.mix.list[[m]] <- model.estimate$latent
      model.observed.list[[m]] <- model.estimate$observed
    }
    
    res.tmp <- intrinsic.variability(y.true.frame = y.true.frame, latent.mix.list = latent.mix.list, 
                                     model.observed.list = model.observed.list, n.samp = 5, N = N,
                                     ker = ker.tmp, h = h.tmp, parallel = T, show.plot = F)
    
    return(res.tmp)
  })
  res.set <- unlist(res.list)
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    l = hyper.param.idx[q,3]
    results.array[i,j,k,l] <- res.set[q]
  }
  
}


#saveRDS(results.array,paste0("data/intrinsic_variability_sims_100.rds"))
results.array <- readRDS("data/intrinsic_variability_sims_100.rds")









# Appendix Simulations 

# true model parameters
ker.model <- ker.set[[1]] #Gaussian Kernel 
h.model <- h.set[7] # 2.0 



beta.par.1 <- 12 
beta.par.2 <- 5
dataset.size <- 100









