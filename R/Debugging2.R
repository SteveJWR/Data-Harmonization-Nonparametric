







for(j in 1:2){
  cond.tmp <- cond.list[[j]]
  A.mat.tmp <- compute_A_matrix(R_bins = R_bins, cond = cond.tmp)
  # replacing error with some kernels
  A.mat.tmp[is.na(A.mat.tmp)] <- 0
  ce.samp.mu.vec <- sapply(1:2, function(l){
    mu.tmp <- mu.set.conversion[l]
    
    # here we introduce a shortcut for computing. We only used the unique elements of the list
    # so we don't have to re-estimate mixing distributions that will be the same. 
    p.hat.list.short <- unique(p.hat.list.z.cw)
    latent.mix.list.short <- list()
    for(k in 1:length(p.hat.list.short)){
      latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.tmp, mu = mu.tmp)
      latent.mix.list.short[[k]] <- latent.model$latent
    }
    
    latent.mix.z.list.cw <- list()
    for(i in 1:length(p.hat.list.z.cw)){
      p.hat.tmp <- p.hat.list.z.cw[[i]]
      short.idx <- which(sapply(p.hat.list.short, function(s) all(s == p.hat.tmp)))
      latent.mix.z.list.cw[[i]] <- latent.mix.list.short[[short.idx]]
    }
    
    
    # sample cross entropy 
    ce.samp  <- convert_score_metric(test.pairs = cw.test, 
                                     latent.mix.list.y = latent.mix.y.list.cw, 
                                     latent.mix.list.z = latent.mix.z.list.cw, 
                                     cond.y = cond.y.opt, 
                                     cond.z = cond.tmp, # using the current iteration of cond.z 
                                     grid.size = 1000) 
    
    
  })
  conversion.results.array[j,1:2] <- ce.samp.mu.vec
  
  ### ML (unregularized estimate)
  
  p.hat.list.short <- unique(p.hat.list.z.cw)
  latent.mix.list.short.ml <- list()
  for(k in 1:length(p.hat.list.short)){
    latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.tmp, mu = 0)
    latent.mix.list.short.ml[[k]] <- latent.model$latent
  }
  
  latent.mix.z.list.cw.ml <- list()
  for(i in 1:length(p.hat.list.z.cw)){
    p.hat.tmp <- p.hat.list.z.cw[[i]]
    short.idx <- which(sapply(p.hat.list.short, function(s) all(s == p.hat.tmp)))
    latent.mix.z.list.cw.ml[[i]] <- latent.mix.list.short.ml[[short.idx]]
  }
  
  
  # sample cross entropy 
  ce.samp.ml <-  convert_score_metric(test.pairs = cw.test, 
                                      latent.mix.list.y = latent.mix.y.list.cw.ml, 
                                      latent.mix.list.z = latent.mix.z.list.cw.ml, 
                                      cond.y = cond.y.opt, 
                                      cond.z = cond.tmp, # using the current iteration of cond.z 
                                      grid.size = 1000) # TODO: Replace with larger value
  
  conversion.results.array.ml[j] <- ce.samp.ml
  
  
  ## Parametric Model 
  z.param <- logitnorm_model_fit(z.train.wide[,1], X = z.train.wide[,2:6], cond = cond.tmp)
  
  ce.param <- compute_conversion_parametric_ce(test.pairs = cw.test.wide[,c(1,2)], X = cw.test.wide[,3:7],
                                               params.y = y.param, params.z = z.param, 
                                               cond.y = cond.y.opt, cond.z = cond.tmp, 
                                               grid.size = 1000)
  
  conversion.results.array.para[j] <- ce.param
  
  cat(paste0("Model: ", j, "/", length(cond.list), " complete" ), end = "\r")
  
}
