


convert_score_metric(test.grid, 
                     latent.mix.list.y, 
                     latent.mix.list.z, 
                     cond.y, 
                     cond.z, 
                     joint.prob = p.yz,
                     grid.size = 1000) 

test.pairs = test.grid
latent.mixture.y = latent.mix.list.y
latent.mixture.z = latent.mix.list.z
cond.y = cond.model1 
cond.z = cond.tmp 
joint.prob = p.yz 
grid.size = 1000







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
  
  pz.given.y <- p.yz[y + 1,]/sum(p.yz[y + 1,])
  
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






