library(parallel)
numCores <- detectCores()
numCores
source("R/01_functions.R")
source("R/02a_simulations_setup.R")


# setting the seed across simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1)






RUN_PARALLEL = FALSE
timeout = 5000
n.sims <- 3




##### ----- intrinsic variability simulations model 1 -----

intrinsic.variability.model1.simulation.results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))



if (RUN_PARALLEL) {
  # Detect number of cores, use all but 1
  no_cores <- detectCores() - 1
  # Initiate cluster
  tictoc::tic()
  
  cl <- makeCluster(no_cores, type="FORK",timeout=timeout)
  
  
  # Run computation
  intrinsic.variability.model1.simulation.results = parLapply(cl = cl, X = 1:n.sims,
                     fun = simulation_intrinsic_variability_model1)
  # Stop cluster
  stopCluster(cl)
  tictoc::toc()
} else {
  tictoc::tic()
  intrinsic.variability.model1.simulation.results = lapply(X = 1:n.sims,
                  FUN = simulation_intrinsic_variability_model1)
  tictoc::toc()
}


for(i in 1:n.sims){
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    l = hyper.param.idx[q,3]
    intrinsic.variability.model1.simulation.results.array[i,j,k,l] <- intrinsic.variability.model1.simulation.results[[i]][q]
  }
}



saveRDS(intrinsic.variability.model1.simulation.results.array,paste0("Data/results/intrinsic_variability_simulation_model1.rds"))

print("intrinsic variability simulations model 1 complete")




##### ----- intrinsic variability simulations model 2 -----

intrinsic.variability.model2.simulation.results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))



if (RUN_PARALLEL) {
  # Detect number of cores, use all but 1
  no_cores <- detectCores() - 1
  # Initiate cluster
  tictoc::tic()
  
  cl <- makeCluster(no_cores, type="FORK",timeout=timeout)
  
  
  # Run computation
  intrinsic.variability.model2.simulation.results = parLapply(cl = cl, X = 1:n.sims,
                                                              fun = simulation_intrinsic_variability_model2)
  # Stop cluster
  stopCluster(cl)
  tictoc::toc()
} else {
  tictoc::tic()
  intrinsic.variability.model2.simulation.results = lapply(X = 1:n.sims,
                                                           FUN = simulation_intrinsic_variability_model2)
  tictoc::toc()
}


for(i in 1:n.sims){
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    l = hyper.param.idx[q,3]
    intrinsic.variability.model2.simulation.results.array[i,j,k,l] <- intrinsic.variability.model2.simulation.results[[i]][q]
  }
}



saveRDS(intrinsic.variability.model2.simulation.results.array,paste0("Data/results/intrinsic_variability_simulation_model2.rds"))

print("intrinsic variability simulations model 2 complete")





## -------- mu selection --------
## --------    Model 1   --------



smoothing.selection.model1.simulation.results.array <- array(NA, dim = c(n.sims, length(mu.set.long)))


if (RUN_PARALLEL) {
  # Detect number of cores, use all but 1
  no_cores <- detectCores() - 1
  # Initiate cluster
  tictoc::tic()
  
  cl <- makeCluster(no_cores, type="FORK",timeout=timeout)
  
  
  # Run computation
  smoothing.selection.model1.simulation.results = parLapply(cl = cl, X = 1:n.sims,
                                                              fun = simulation_smoothing_selection_model1)
  # Stop cluster
  stopCluster(cl)
  tictoc::toc()
} else {
  tictoc::tic()
  smoothing.selection.model1.simulation.results = lapply(X = 1:n.sims,
                                                           FUN = simulation_smoothing_selection_model1)
  tictoc::toc()
}


for(i in 1:n.sims){
  for(q in 1:length(mu.set.long)){
    smoothing.selection.model1.simulation.results.array[i,q] <- smoothing.selection.model1.simulation.results[[i]][q]
  }
}



saveRDS(smoothing.selection.model1.simulation.results.array,paste0("Data/results/smoothing_selection_simulation_model1.rds"))


print("mu selection simulations model 1 complete")


## --------    Model 2   --------


smoothing.selection.model2.simulation.results.array <- array(NA, dim = c(n.sims, length(mu.set.long)))


if (RUN_PARALLEL) {
  # Detect number of cores, use all but 1
  no_cores <- detectCores() - 1
  # Initiate cluster
  tictoc::tic()
  
  cl <- makeCluster(no_cores, type="FORK",timeout=timeout)
  
  # Run computation
  smoothing.selection.model2.simulation.results = parLapply(cl = cl, X = 1:n.sims,
                                                            fun = simulation_smoothing_selection_model2)
  # Stop cluster
  stopCluster(cl)
  tictoc::toc()
} else {
  tictoc::tic()
  smoothing.selection.model2.simulation.results = lapply(X = 1:n.sims,
                                                         FUN = simulation_smoothing_selection_model2)
  tictoc::toc()
}


for(i in 1:n.sims){
  for(q in 1:length(mu.set.long)){
    smoothing.selection.model2.simulation.results.array[i,q] <- smoothing.selection.model2.simulation.results[[i]][q]
  }
}



saveRDS(smoothing.selection.model2.simulation.results.array,paste0("Data/results/smoothing_selection_simulation_model2.rds"))


print("mu selection simulations model 2 complete")





## -------- Conversion Cross Entropy --------



conversion.ce.simulation.results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set.conversion)))
conversion.ce.simulation.ml.results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set)))



if (RUN_PARALLEL) {
  # Detect number of cores, use all but 1
  no_cores <- detectCores() - 1
  # Initiate cluster
  tictoc::tic()
  
  cl <- makeCluster(no_cores, type="FORK",timeout=timeout)
  
  
  # Run computation
  conversion.ce.simulation.results = parLapply(cl = cl, X = 1:n.sims,
                                               fun = simulation_conversion_cross_entropy)
  # Stop cluster
  stopCluster(cl)
  tictoc::toc()
} else {
  tictoc::tic()
  conversion.ce.simulation.results = lapply(X = 1:n.sims,
                                            FUN = simulation_conversion_cross_entropy)
  tictoc::toc()
}


for(i in 1:n.sims){
  for(q in 1:nrow(hyper.param.conversion.idx)){
    j = hyper.param.conversion.idx[q,1]
    k = hyper.param.conversion.idx[q,2]
    l = hyper.param.conversion.idx[q,3]
    conversion.ce.simulation.results.array[i,j,k,l] <- conversion.ce.simulation.results[[i]]$smoothed[q]
  }
  for(q in 1:nrow(hyper.param.conversion.ml.idx)){
    j = hyper.param.conversion.ml.idx[q,1]
    k = hyper.param.conversion.ml.idx[q,2]
    conversion.ce.simulation.ml.results.array[i,j,k] <- conversion.ce.simulation.results[[i]]$npmle[q]
  }
}

saveRDS(conversion.ce.simulation.results.array,paste0("Data/results/conversion_ce_simulations.rds"))
saveRDS(conversion.ce.simulation.ml.results.array,paste0("Data/results/conversion_ce_ml_simulations.rds"))

print("conversion cross entropy simulations complete")






## -------- Feasibility Tests -------- 














