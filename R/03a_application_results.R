
# TODO: DELETE 
# setwd("/Users/Owner/Documents/PhD/Latent Variable Modelling with Application to Data Harmonization/GitHub/Data-Harmonization-Nonparametric")
source("R/01_functions.R")

# TODO: Add these to the list of functions
library(dplyr)
library(gridExtra)
library(ggpubr)

####
##### Figure information 
png.width = 1200
png.height = 1000

title.size = 40 #22
axis.ticks.size = 24 #14
axis.size = 30 #18

####


y.train <- read.csv("Data/NACCMMSE_training.csv")
z.train <- read.csv("Data/MOCATOTS_training.csv")

y.val <- read.csv("Data/NACCMMSE_validation.csv")
z.val <- read.csv("Data/MOCATOTS_validation.csv")

y.val <- y.val[,c("y1", "y2", "age", "group")]
z.val <- z.val[,c("z1", "z2", "age", "group")]


Ny <- 30
Nz <- 30




h.set<- exp(seq(log(0.5),log(15), length.out = 15))
ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
ker.names <- c("Gaussian","Exponential", "Triangle", "Epanechnikov")
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 2)))

### Perhaps alter this section for the real data analysis 
results.array <- array(NA, dim = c(length(h.set), length(ker.set), length(mu.set)))
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))


# remove unregularized compact support pairs (causes errors)
idx1 <- hyper.param.idx[,1] <= 1.0  
idx2 <- hyper.param.idx[,2] %in% c(3,4) 

# Also remove compact kernels with h < 1
idx3 <- hyper.param.idx[,1] < 1 

hyper.param.idx <- hyper.param.idx[!(idx1 & idx2 ), ]



## Transforming the data set to the corresponding conditional estimators 



age.range <- 60:85
group.range <- 1:4
x.design <- data.frame(age = rep(age.range, times = length(group.range)), 
                       group = rep(group.range, each = length(age.range)))



# computing a list of conditional estimators corresponding to the design matrix. 

p.hat.short.list.y <- list()
for(j in 1:nrow(x.design)){
  age.tmp <- x.design$age[j]
  group.tmp <- x.design$group[j]
  p.hat.tmp <- estimate_conditional_distribution(train.data = y.train, 
                                                 group = group.tmp, 
                                                 age = age.tmp, 
                                                 outcome = "y", 
                                                 N = Ny)
  
  #### 
  p.hat.short.list.y[[j]] <-   p.hat.tmp
}



p.hat.short.list.z <- list()
for(j in 1:nrow(x.design)){
  age.tmp <- x.design$age[j]
  group.tmp <- x.design$group[j]
  p.hat.tmp <- estimate_conditional_distribution(train.data = z.train, 
                                                 group = group.tmp, 
                                                 age = age.tmp, 
                                                 outcome = "z", 
                                                 N = Nz)
  
  #### 
  p.hat.short.list.z[[j]] <-   p.hat.tmp
}



####### Additional Functions 
estimate_conditional_distribution <- function(train.data, group, age, outcome, N, age.window = 3){
  age.idx <- which(colnames(train.data) == "age")
  group.idx <- which(colnames(train.data) == "group")
  filter.idx <- abs(as.vector(train.data[,age.idx]) - age) <= 3 & (as.vector(train.data[,group.idx]) == group)
  data.subset <- train.data[filter.idx,]
  p.hat <- compute_edf(data.subset[,outcome], N)
  
  return(p.hat)
}


# 
sort_conditionals <- function(data, x.design, p.hat.list){
  n = nrow(data)
  out <- list()
  for(i in 1:n){
    age.tmp <- data$age[j]
    group.tmp <- data$group[j]
    
    j <- which( x.design$age == age.tmp & x.design$group == group.tmp)
    
    out[[i]] <- p.hat.list[[j]]
  }
  return(out)
}







# include a function which already matches the p.hat to each row, pair.obs only includes outcomes. 

select_cond <- function(pair.obs, p.hat, cond.list, cond.names, mu, n.samp = 10, R_bins = 1000){
  n <- nrow(pair.obs)
  
  if(missing(cond.names) | length(cond.list) != length(cond.names)){
    cond.names <- as.character(1:length(cond.list))
  }
  
  # a resonable default found in practice 
  if(missing(mu)){
    mu = 0.1/(length(cond(0.5)))
  }
  
  
  # check if the current p.hat is a "list"
  if(typeof(p.hat) == "numeric" ){
    p.hat.list <- list()
    for(i in 1:nrow(pair.obs)){
      p.hat.list[[i]] <- p.hat
    }
  } else if(nrow(pair.obs) == length(p.hat) & typeof(p.hat) == "list"){
    p.hat.list <- p.hat
  } else {
    stop("p.hat must be either a vector or a list of the same length as the observations")
  }
  
  # vector of values for the intrinsic variability
  iv.dists <- rep(NA, length(cond.list))
  for(j in 1:length(cond.list)){
    cond.tmp <- cond.list[[j]]
    A.mat.tmp <- compute_A_matrix(R_bins = R_bins, cond = cond.tmp)
    
    
    # here we introduce a shortcut for computing. We only used the unique elements of the list
    # so we don't have to re-estimate mixing distributions that will be the same. 
    p.hat.list.short <- unique(p.hat.list)
    latent.mix.list.short <- list()
    for(k in 1:length(p.hat.list.short)){
      latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.tmp, mu = mu)
      latent.mix.list.short[[k]] <- latent.model$latent
    }
    
    latent.mix.list <- list()
    for(i in 1:n){
      p.hat.tmp <- p.hat.list[[i]]
      short.idx <- which(sapply(p.hat.list.short, function(s) all(match(s, p.hat.tmp))))
      latent.mix.list[[i]] <- latent.mix.list.short[[short.idx]]
    }
    iv.dists[j] <- intrinsic_variability(pair.obs, latent.mix.list, n.samp, cond.tmp)
    
    cat(paste0("Model: ", j, "/", length(cond.list), " complete" ), end = "\r")
  }
  opt.idx <- which.min(iv.dists)
  
  print(paste0("Optimal Conditional Model: ", cond.names[opt.idx], "|| TV dist: ", min(iv.dists)))
  out.list <- list("results" = iv.dists, "optimal" = cond.list[[opt.idx]])
  return(out.list)
}



select_cond <- function(pair.obs, p.hat, cond.list, cond.names, mu, n.samp = 10, R_bins = 1000){
  n <- nrow(pair.obs)
  
  if(missing(cond.names) | length(cond.list) != length(cond.names)){
    cond.names <- as.character(1:length(cond.list))
  }
  
  # a resonable default found in practice 
  if(missing(mu)){
    mu = 0.1/(length(cond(0.5)))
  }
  
  
  # check if the current p.hat is a "list"
  if(typeof(p.hat) == "numeric" ){
    p.hat.list <- list()
    for(i in 1:nrow(pair.obs)){
      p.hat.list[[i]] <- p.hat
    }
  } else if(nrow(pair.obs) == length(p.hat) & typeof(p.hat) == "list"){
    p.hat.list <- p.hat
  } else {
    stop("p.hat must be either a vector or a list of the same length as the observations")
  }
  
  # vector of values for the intrinsic variability
  iv.dists <- rep(NA, length(cond.list))
  for(j in 1:length(cond.list)){
    cond.tmp <- cond.list[[j]]
    A.mat.tmp <- compute_A_matrix(R_bins = R_bins, cond = cond.tmp)
    
    
    # here we introduce a shortcut for computing. We only used the unique elements of the list
    # so we don't have to re-estimate mixing distributions that will be the same. 
    p.hat.list.short <- unique(p.hat.list)
    latent.mix.list.short <- list()
    for(k in 1:length(p.hat.list.short)){
      latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat.tmp, mu = mu)
      latent.mix.list.short[[k]] <- latent.model$latent
    }
    
    latent.mix.list <- list()
    for(i in 1:n){
      p.hat.tmp <- p.hat.list[[i]]
      short.idx <- which(sapply(p.hat.list.short, function(s) all(match(s, p.hat.tmp))))
      latent.mix.list[[i]] <- latent.mix.list.short[[short.idx]]
    }
    iv.dists[j] <- intrinsic_variability(pair.obs, latent.mix.list, n.samp, cond.tmp)
    
    cat(paste0("Model: ", j, "/", length(cond.list), " complete" ), end = "\r")
  }
  opt.idx <- which.min(iv.dists)
  
  print(paste0("Optimal Conditional Model: ", cond.names[opt.idx], "|| TV dist: ", min(iv.dists)))
  out.list <- list("results" = iv.dists, "optimal" = cond.list[[opt.idx]])
  return(out.list)
}


select_mu <- function(pair.obs, p.hat, cond, mu.vec, R_bins = 1000){
  n <- nrow(pair.obs)
  A.mat <- compute_A_matrix(R_bins = R_bins, cond = cond.tmp)
  A.tensor <- compute_A_two_obs_tensor(R_bins = R_bins, cond = cond.tmp)
  
  
  # check if the current p.hat is a "list"
  if(typeof(p.hat) == "numeric" ){
    p.hat.list <- list()
    for(i in 1:nrow(pair.obs)){
      p.hat.list[[i]] <- p.hat
    }
  } else if(nrow(pair.obs) == length(p.hat) & typeof(p.hat) == "list"){
    p.hat.list <- p.hat
  } else {
    stop("p.hat must be either a vector or a list of the same length as the observations")
  }
  
  
  # vector of values for the two sample likelihood
  lik.vals <- rep(NA, length(cond.list))
  for(j in 1:length(mu.vec)){
    mu.tmp <- mu.vec[j]
    
   
    
    # here we introduce a shortcut for computing. We only used the unique elements of the list
    # so we don't have to re-estimate mixing distributions that will be the same. 
    p.hat.list.short <- unique(p.hat.list)
    latent.mix.list.short <- list()
    for(k in 1:length(p.hat.list.short)){
      latent.model <- estimate_mixing_numeric(p.hat.list.short[[k]], A.mat, mu = mu.tmp)
      latent.mix.list.short[[k]] <- latent.model$latent
    }
    
    latent.mix.list <- list()
    for(i in 1:n){
      p.hat.tmp <- p.hat.list[[i]]
      short.idx <- which(sapply(latent.mix.list.short, function(s) all(match(s, p.hat.tmp))))
      latent.mix.list[[i]] <- latent.mix.list.short[[short.idx]]
    }
    
    lik.vals[j] <- compute_two_obs_loglikelihood(pair.obs, A.tensor, latent.mix.list)
    
    cat(paste0("Smoothing Parameter: ", j, "/", length(mu.vec), " complete" ) , end = "\r")
  }
  opt.idx <- which.max(lik.vals)
  
  print(paste0("Optimal Smoothing Parameter: ", mu.vec[opt.idx]))
  out.list <- list("results" = mu.vec, "optimal" = mu.vec[opt.idx])
  return(out.list)
}


generate_mkm_list <- function(N,ker,h.set){
  out.list <- list()
  for(i in 1:length(h.set)){
    out.list[[i]] <- conditional_mkm(N,ker,h.set[i])
  }
  return(out.list)
}

generate_cond_binomial <- function(N){
  out.function <- function(x){
    dbinom(0:N,size = N, prob = x)
  }
}

##### Model selection: 

gaussian.cond.list <- generate_mkm_list(N = Ny, ker = gaussian_kernel, h.set = h.set)
gaussian.cond.names <- paste0("Gaussian Kernel, h = ", round(h.set,2))

exponential.cond.list <- generate_mkm_list(N = Ny, ker = exponential_kernel, h.set = h.set)
exponential.cond.names <- paste0("Exponential Kernel, h = ", round(h.set,2))

triangle.cond.list <- generate_mkm_list(N = Ny, ker = triangle_kernel, h.set = h.set)
triangle.cond.names <- paste0("Triangle Kernel, h = ", round(h.set,2))

epanechnikov.cond.list <- generate_mkm_list(N = Ny, ker = epanechnikov_kernel, h.set = h.set)
epanechnikov.cond.names <- paste0("Epanechnikov Kernel, h = ", round(h.set,2))


cond.list <- c(gaussian.cond.list,exponential.cond.list,triangle.cond.list,epanechnikov.cond.list, generate_cond_binomial(Ny))
cond.names <- c(gaussian.cond.names,exponential.cond.names,triangle.cond.names,epanechnikov.cond.names, "Binomial")





#### MMSE Model Selection 
p.hat.list.y.val <- sort_conditionals(data = y.val, 
                                      x.design = x.design, 
                                      p.hat.list = p.hat.short.list.y)


### Applying this to three choices of mu for completeness
cond.y.model.0 <- select_cond(pair.obs = y.val,  
                              p.hat = p.hat.list.y.val, 
                              cond.list = cond.list, 
                              cond.names = cond.names, 
                              mu = 0, 
                              n.samp = 10, 
                              R_bins = 1000)

cond.y.model.0.001 <- select_cond(pair.obs = y.val,  
                                 p.hat = p.hat.list.y.val, 
                                 cond.list = cond.list, 
                                 cond.names = cond.names, 
                                 mu = 0.001, 
                                 n.samp = 10, 
                                 R_bins = 1000)

cond.y.model.0.1 <- select_cond(pair.obs = y.val,  
                                p.hat = p.hat.list.y.val, 
                                cond.list = cond.list, 
                                cond.names = cond.names, 
                                mu = 0.1, 
                                n.samp = 10, 
                                R_bins = 1000)



saveRDS(cond.y.model.0, "conditional_model_selection_mu0_mmse.rds")
saveRDS(cond.y.model.0.001, "conditional_model_selection_mu0001_mmse.rds")
saveRDS(cond.y.model.0.1, "conditional_model_selection_mu01_mmse.rds")


#### 



mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(ker.set))

gauss.tv <- c(cond.y.model.0$results[1:length(h.set)],
              cond.y.model.0.001$results[1:length(h.set)],
              cond.y.model.0.1$results[1:length(h.set)])

exp.tv <- c(cond.y.model.0$results[(length(h.set) + 1:length(h.set))],
            cond.y.model.0.001$results[(length(h.set) + 1:length(h.set))],
            cond.y.model.0.1$results[(length(h.set) + 1:length(h.set))])

triangle.tv <- c(cond.y.model.0$results[(2*length(h.set) + 1:length(h.set))],
                 cond.y.model.0.001$results[(2*length(h.set) + 1:length(h.set))],
                 cond.y.model.0.1$results[(2*length(h.set) + 1:length(h.set))])

epanechnikov.tv <- c(cond.y.model.0$results[(3*length(h.set) + 1:length(h.set))],
                     cond.y.model.0.001$results[(3*length(h.set) + 1:length(h.set))],
                     cond.y.model.0.1$results[(3*length(h.set) + 1:length(h.set))])


gauss.frame <- data.frame(h = h.set, mu = mus, total.variation = gauss.tv)
exp.frame <- data.frame(h = h.set, mu = mus, total.variation = exp.tv)
triangle.frame <- data.frame(h = h.set, mu = mus, total.variation = triangle.tv)
epanechnikov.frame <- data.frame(h = h.set, mu = mus, total.variation = epanechnikov.tv)


binom.tv <- cond.y.model.0.001$results[(4*length(h.set) + 1)]

# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Gaussian Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                        axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                        axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                        axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                        title = element_text( size = title.size),
                                                                                                                                                                                                                        legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                        legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                        legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                        legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                        legend.text = element_text(size=axis.ticks.size))
exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Laplace Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                   axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                   axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                   axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                   title = element_text( size = title.size),
                                                                                                                                                                                                                   legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                   legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                   legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                   legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                   legend.text = element_text(size=axis.ticks.size))
triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Triangle Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                               axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                               axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                               axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                               title = element_text( size = title.size),
                                                                                                                                                                                                                               legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                               legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                               legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                               legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                               legend.text = element_text(size=axis.ticks.size))
epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Epanechnikov Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                           axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                           axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                           axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                           title = element_text( size = title.size),
                                                                                                                                                                                                                                           legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                           legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                           legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                           legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                           legend.text = element_text(size=axis.ticks.size))



png(filename = "Plots/Intrinsic_Variability_MMSE.png",
    width = png.width, height = png.height)

ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 





##### MOCA Model Selection 


p.hat.list.z.val <- sort_conditionals(data = z.val, 
                                      x.design = x.design, 
                                      p.hat.list = p.hat.short.list.z)


### Applying this to three choices of mu for completeness
cond.z.model.0 <- select_cond(pair.obs = z.val,  
                              p.hat = p.hat.list.z.val, 
                              cond.list = cond.list, 
                              cond.names = cond.names, 
                              mu = 0, 
                              n.samp = 10, 
                              R_bins = 1000)

cond.z.model.0.001 <- select_cond(pair.obs = z.val,  
                                  p.hat = p.hat.list.z.val, 
                                  cond.list = cond.list, 
                                  cond.names = cond.names, 
                                  mu = 0.001, 
                                  n.samp = 10, 
                                  R_bins = 1000)

cond.z.model.0.1 <- select_cond(pair.obs = z.val,  
                                p.hat = p.hat.list.z.val, 
                                cond.list = cond.list, 
                                cond.names = cond.names, 
                                mu = 0.1, 
                                n.samp = 10, 
                                R_bins = 1000)



saveRDS(cond.z.model.0, "conditional_model_selection_mu0_moca.rds")
saveRDS(cond.z.model.0.001, "conditional_model_selection_mu0001_moca.rds")
saveRDS(cond.z.model.0.1, "conditional_model_selection_mu01_moca.rds")


#### 



mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(ker.set))

gauss.tv <- c(cond.z.model.0$results[1:length(h.set)],
              cond.z.model.0.001$results[1:length(h.set)],
              cond.z.model.0.1$results[1:length(h.set)])

exp.tv <- c(cond.z.model.0$results[(length(h.set) + 1:length(h.set))],
            cond.z.model.0.001$results[(length(h.set) + 1:length(h.set))],
            cond.z.model.0.1$results[(length(h.set) + 1:length(h.set))])

triangle.tv <- c(cond.z.model.0$results[(2*length(h.set) + 1:length(h.set))],
                 cond.z.model.0.001$results[(2*length(h.set) + 1:length(h.set))],
                 cond.z.model.0.1$results[(2*length(h.set) + 1:length(h.set))])

epanechnikov.tv <- c(cond.z.model.0$results[(3*length(h.set) + 1:length(h.set))],
                     cond.z.model.0.001$results[(3*length(h.set) + 1:length(h.set))],
                     cond.z.model.0.1$results[(3*length(h.set) + 1:length(h.set))])


gauss.frame <- data.frame(h = h.set, mu = mus, total.variation = gauss.tv)
exp.frame <- data.frame(h = h.set, mu = mus, total.variation = exp.tv)
triangle.frame <- data.frame(h = h.set, mu = mus, total.variation = triangle.tv)
epanechnikov.frame <- data.frame(h = h.set, mu = mus, total.variation = epanechnikov.tv)


binom.tv <- cond.z.model.0.001$results[(4*length(h.set) + 1)]

# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Gaussian Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                        axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                        axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                        axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                        title = element_text( size = title.size),
                                                                                                                                                                                                                        legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                        legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                        legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                        legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                        legend.text = element_text(size=axis.ticks.size))
exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Laplace Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                   axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                   axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                   axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                   title = element_text( size = title.size),
                                                                                                                                                                                                                   legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                   legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                   legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                   legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                   legend.text = element_text(size=axis.ticks.size))
triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Triangle Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                               axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                               axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                               axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                               title = element_text( size = title.size),
                                                                                                                                                                                                                               legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                               legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                               legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                               legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                               legend.text = element_text(size=axis.ticks.size))
epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Epanechnikov Kernel") + geom_hline(yintercept = binom.tv) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                           axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                           axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                           axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                           title = element_text( size = title.size),
                                                                                                                                                                                                                                           legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                           legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                           legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                           legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                           legend.text = element_text(size=axis.ticks.size))



png(filename = "Plots/Intrinsic_Variability_MOCA.png",
    width = png.width, height = png.height)

ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 
















## the rest of the analysis will just be 
## list the grid of conditional distributions 
## plot the results for the conditional and mu selection procedures. 
## Then plot the conversion probability as a function of the second model 
## Also must include the conversion cross entropy for the parametric as well as simple z-score conversion methods 




