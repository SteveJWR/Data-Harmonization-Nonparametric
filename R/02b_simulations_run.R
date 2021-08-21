library(parallel)
numCores <- detectCores()
numCores
source("R/01_functions.R")
source("R/02a_simulations_setup.R")


# setting the seed across simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1)



##### ----- intrinsic variability simulations model 1 -----

intrinsic.variability.model1.simulation.results.array <- array(NA, dim = c(n.sims, length(h.set), length(ker.set), length(mu.set)))
intrinsic.variability.model1.simulation.results <- mclapply(1:n.sims,
                                                            simulation_intrinsic_variability_model1)


for(i in 1:n.sims){
  for(q in 1:nrow(hyper.param.idx1)){
    j = hyper.param.idx1[q,1]
    k = hyper.param.idx1[q,2]
    l = hyper.param.idx1[q,3]
    intrinsic.variability.model1.simulation.results.array[i,j,k,l] <- intrinsic.variability.model1.simulation.results[[i]][q]
  }
}



saveRDS(intrinsic.variability.model1.simulation.results.array,paste0("data/intrinsic_variability_simulation_model1.rds"))





