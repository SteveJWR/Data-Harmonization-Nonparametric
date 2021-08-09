

setwd("/Users/Owner/Documents/PhD/Latent Variable Modelling with Application to Data Harmonization")

# Functions for Data Harmonization
source("DataHarmonizationFunctions.R")
# Note that we must install the MOSEK solver as well as the MultChernoff bound 

# require this in the packages 
library(gridExtra)
library(ggpubr)



##### Figure information 
png.width = 1200
png.height = 1000

title.size = 40 #22
axis.ticks.size = 24#14
axis.size = 30 #18






h.set<- c(0.2,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5, 8.0) 

ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 3)))


# true model parameters
ker.model <- ker.set[[1]] #Gaussian Kernel 
h.model <- h.set[7] # 2.0 



beta.par.1 <- 12 
beta.par.2 <- 5
dataset.size <- 100





set.seed(2)
N = 30
dataset.size = 100000
A.model <- A.matrix.compute(R_bins = 1000, N = 30, ker = ker.model, h = h.model, numeric.points = 400)
latent.model.samp <- rbeta(dataset.size, shape1 = beta.par.1, shape2 = beta.par.2)

y.model.samp <- lapply(latent.model.samp, function(x){
  i = 0:N
  model.weights <- ker.model((i - N*x)/h.model)
  out <- sample(0:N, size = 1, replace = TRUE, prob = model.weights)
  
})


y.model.samp <- unlist(y.model.samp)

true.p.hat <- c()
for(y in 0:(N)){
  prop.tmp <- mean(y.model.samp == y)
  true.p.hat <- c(true.p.hat, prop.tmp)
}



set.seed(2)
N = 30
dataset.size = 100
R_bins = 1000
A.model <- A.matrix.compute(R_bins = R_bins, N = 30, ker = ker.model, h = h.model, numeric.points = 400)
latent.model.samp <- rbeta(dataset.size, shape1 = beta.par.1, shape2 = beta.par.2)

y.model.samp <- lapply(latent.model.samp, function(x){
  i = 0:N
  model.weights <- ker.model((i - N*x)/h.model)
  out <- sample(0:N, size = 1, replace = TRUE, prob = model.weights)
  
})


y.model.samp <- unlist(y.model.samp)

train.p.hat <- c()
for(y in 0:(N)){
  prop.tmp <- mean(y.model.samp == y)
  train.p.hat <- c(train.p.hat, prop.tmp)
}






# Visualising different mu values 
########################
mu = 0.1 # 0.001, 0.01, 0.1, 0.4




model.estimate0 <- numeric.latent.fit(p.hat = train.p.hat, 
                                      A.matrix = A.model, 
                                      mu = 0, show.plot = F)

model.estimate0.001 <- numeric.latent.fit(p.hat = train.p.hat, 
                                          A.matrix = A.model, 
                                          mu = 0.001, show.plot = F)

model.estimate0.01 <- numeric.latent.fit(p.hat = train.p.hat, 
                                         A.matrix = A.model, 
                                         mu = 0.01, show.plot = F)

model.estimate0.1 <- numeric.latent.fit(p.hat = train.p.hat, 
                                        A.matrix = A.model, 
                                        mu = 0.1, show.plot = F) 

model.estimate0.4 <- numeric.latent.fit(p.hat = train.p.hat, 
                                        A.matrix = A.model, 
                                        mu = 0.4, show.plot = F) 

gamma.seq <- seq(0,1,length.out = R_bins)

cdf.0 <- rep(NA, R_bins)
cdf.0.001 <- rep(NA, R_bins)
cdf.0.01 <- rep(NA, R_bins)
cdf.0.1 <- rep(NA, R_bins)
cdf.0.4 <- rep(NA, R_bins)

for(k in 1:R_bins){
  cdf.0[k] <- sum(model.estimate0$latent[1:k])
  cdf.0.001[k] <- sum(model.estimate0.001$latent[1:k])
  cdf.0.01[k] <- sum(model.estimate0.01$latent[1:k])
  cdf.0.1[k] <- sum(model.estimate0.1$latent[1:k])
  cdf.0.4[k] <- sum(model.estimate0.4$latent[1:k])
}


latent.true <- dbeta(gamma.seq, beta.par.1, beta.par.2)
cdf.true <- pbeta(gamma.seq, beta.par.1, beta.par.2)

y.seq = 0:N



### Observed Data
obs.mu0 <- data.frame(y = rep(y.seq, 3), 
                      labels = rep(c("True", "Empirical", "Estimate"), 
                                   each = length(y.seq)),
                      probability = c(true.p.hat, 
                                      train.p.hat, 
                                      as.vector(model.estimate0$observed)))


obs.mu0.001 <- data.frame(y = rep(y.seq, 3), 
                      labels = rep(c("True", "Empirical", "Estimate"), 
                                   each = length(y.seq)),
                      probability = c(true.p.hat, 
                                      train.p.hat, 
                                      as.vector(model.estimate0.001$observed)))


obs.mu0.01 <- data.frame(y = rep(y.seq, 3), 
                      labels = rep(c("True", "Empirical", "Estimate"), 
                                   each = length(y.seq)),
                      probability = c(true.p.hat, 
                                      train.p.hat, 
                                      as.vector(model.estimate0.01$observed)))


obs.mu0.1 <- data.frame(y = rep(y.seq, 3), 
                      labels = rep(c("True", "Empirical", "Estimate"), 
                                   each = length(y.seq)),
                      probability = c(true.p.hat, 
                                      train.p.hat, 
                                      as.vector(model.estimate0.1$observed)))

obs.mu0.4 <- data.frame(y = rep(y.seq, 3), 
                        labels = rep(c("True", "Empirical", "Estimate"), 
                                     each = length(y.seq)),
                        probability = c(true.p.hat, 
                                        train.p.hat, 
                                        as.vector(model.estimate0.4$observed)))


### Latent Density 


dens.mu0 <- data.frame(y = rep(gamma.seq, 2), 
                      labels = rep(c("True", "Estimate"), 
                                   each = length(gamma.seq)),
                      density = c(latent.true, 
                                      R_bins*as.vector(model.estimate0$latent)))


dens.mu0.001 <- data.frame(y = rep(gamma.seq, 2), 
                           labels = rep(c("True", "Estimate"), 
                                        each = length(gamma.seq)),
                           density = c(latent.true, 
                                       R_bins*as.vector(model.estimate0.001$latent)))

dens.mu0.01 <- data.frame(y = rep(gamma.seq, 2), 
                          labels = rep(c("True", "Estimate"), 
                                       each = length(gamma.seq)),
                          density = c(latent.true, 
                                      R_bins*as.vector(model.estimate0.01$latent)))

dens.mu0.1 <- data.frame(y = rep(gamma.seq, 2), 
                         labels = rep(c("True", "Estimate"), 
                                      each = length(gamma.seq)),
                         density = c(latent.true, 
                                     R_bins*as.vector(model.estimate0.1$latent)))

dens.mu0.4 <- data.frame(y = rep(gamma.seq, 2), 
                         labels = rep(c("True", "Estimate"), 
                                      each = length(gamma.seq)),
                         density = c(latent.true, 
                                     R_bins*as.vector(model.estimate0.4$latent)))




### CDF 


cdf.mu0 <-  data.frame(y = rep(gamma.seq, 2), 
                       labels = rep(c("True", "Estimate"), 
                                    each = length(gamma.seq)),
                       cdf = c(cdf.true, 
                               cdf.0))

cdf.mu0.001 <- data.frame(y = rep(gamma.seq, 2), 
                          labels = rep(c("True", "Estimate"), 
                                       each = length(gamma.seq)),
                          cdf = c(cdf.true, 
                                  cdf.0.001))

cdf.mu0.01 <- data.frame(y = rep(gamma.seq, 2), 
                         labels = rep(c("True", "Estimate"), 
                                      each = length(gamma.seq)),
                         cdf = c(cdf.true,
                                 cdf.0.01))

cdf.mu0.1 <- data.frame(y = rep(gamma.seq, 2), 
                        labels = rep(c("True", "Estimate"), 
                                     each = length(gamma.seq)),
                        cdf = c(cdf.true,
                                cdf.0.1))

cdf.mu0.4 <- data.frame(y = rep(gamma.seq, 2), 
                        labels = rep(c("True", "Estimate"), 
                                     each = length(gamma.seq)),
                        cdf = c(cdf.true, 
                                cdf.0.4))




p.obs.0 <- ggplot(data = obs.mu0, aes(x = y, 
                                      y = probability, 
                                      group = labels, 
                                      color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 


p.obs.0.001 <- ggplot(data = obs.mu0.001, aes(x = y, 
                                      y = probability, 
                                      group = labels, 
                                      color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.001)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 



p.obs.0.01 <- ggplot(data = obs.mu0.01, aes(x = y, 
                                      y = probability, 
                                      group = labels, 
                                      color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.01)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 


p.obs.0.1 <- ggplot(data = obs.mu0.1, aes(x = y, 
                                      y = probability, 
                                      group = labels, 
                                      color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.1)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 


p.obs.0.4 <- ggplot(data = obs.mu0.4, aes(x = y, 
                                        y = probability, 
                                        group = labels, 
                                        color = labels)) + 
  geom_point() + scale_colour_manual(values = c("red", "blue", "black")) + 
  ggtitle(paste0("Observed Score Estimate \u00b5 = ", 0.4)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 



##### Latent


p.dens.0 <- ggplot(data = dens.mu0, aes(x = y, 
                                      y = density, 
                                      group = labels, 
                                      color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) #+ 
  # labs(col="Density") + annotate(geom="text", x=0.1, y= 500, label=paste0("\u00b5 = ", 0),
  #                                color="blue", size = 6)


p.dens.0


p.dens.0.001 <- ggplot(data = dens.mu0.001, aes(x = y, 
                                          y = density, 
                                          group = labels, 
                                          color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.001)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) #+ 
  # labs(col="Density") + annotate(geom="text", x=0.1, y= 45, label=paste0("\u00b5 = ", 0.001),
  #                                color="blue", size = 6)

p.dens.0.001



p.dens.0.01 <- ggplot(data = dens.mu0.01, aes(x = y, 
                                         y = density, 
                                         group = labels, 
                                         color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.01)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) #+ 
  # labs(col="Density") + annotate(geom="text", x=0.1, y= 12, label=paste0("\u00b5 = ", 0.01),
  #                                color="blue", size = 6)

p.dens.0.01


p.dens.0.1 <- ggplot(data = dens.mu0.1, aes(x = y, 
                                        y = density, 
                                        group = labels, 
                                        color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.1)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) #+ 
  # labs(col="Density") + annotate(geom="text", x=0.1, y= 6, label=paste0("\u00b5 = ", 0.1),
  #                                color="blue", size = 6)

p.dens.0.1

p.dens.0.4 <- ggplot(data = dens.mu0.4, aes(x = y, 
                                           y = density, 
                                           group = labels, 
                                           color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("Latent Density Estimate \u00b5 = ", 0.4)) + 
  xlab("Gamma") + ylab("Density") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) #+ 
  # labs(col="Density") + annotate(geom="text", x=0.1, y= 4, label=paste0("\u00b5 = ", 0.4),
  #            color="blue", size = 6)

p.dens.0.4


##### CDF


p.cdf.0 <- ggplot(data = cdf.mu0, aes(x = y, 
                                        y = cdf, 
                                        group = labels, 
                                        color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 


p.cdf.0.001 <- ggplot(data = cdf.mu0.001, aes(x = y, 
                                               y = cdf, 
                                               group = labels, 
                                               color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.001)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 



p.cdf.0.01 <- ggplot(data = cdf.mu0.01, aes(x = y, 
                                             y = cdf, 
                                             group = labels, 
                                             color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.01)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 


p.cdf.0.1 <- ggplot(data = cdf.mu0.1, aes(x = y, 
                                           y = cdf, 
                                           group = labels, 
                                           color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.1)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size)) 


p.cdf.0.4 <- ggplot(data = cdf.mu0.4, aes(x = y, 
                                          y = cdf, 
                                          group = labels, 
                                          color = labels)) + 
  geom_line() + scale_colour_manual(values = c("blue", "black")) + 
  ggtitle(paste0("CDF Estimate \u00b5 = ", 0.4)) + xlab("Gamma")+ 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  






png(filename = "plots/Example Regularization.png",
    width = (1.5)*png.width, height = png.height)
# 
# grid.arrange(arrangeGrob(p.obs.0, p.dens.0, p.cdf.0, ncol = 3),                             
#              arrangeGrob(p.obs.0.001, p.dens.0.001, p.cdf.0.001, ncol = 3), 
#              arrangeGrob(p.obs.0.01, p.dens.0.01, p.cdf.0.01, ncol = 3), 
#              arrangeGrob(p.obs.0.1, p.dens.0.1, p.cdf.0.1, ncol = 3), 
#              arrangeGrob(p.obs.0.4, p.dens.0.4, p.cdf.0.4, ncol = 3), 
#              nrow = 5)   
ggarrange(p.obs.0, p.dens.0, p.cdf.0,
          p.obs.0.001, p.dens.0.001, p.cdf.0.001,
          p.obs.0.01, p.dens.0.01, p.cdf.0.01,
          p.obs.0.1, p.dens.0.1, p.cdf.0.1,
          p.obs.0.4, p.dens.0.4, p.cdf.0.4,
          ncol=3, nrow=5, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 



####### Plots for AAIC
p.dens.0
p.dens.0.01
p.dens.0.1
p.dens.0.4









######################



############### Speed Test ##############


######## SIMULATIONS ######### 



n.sample.sims = 100
dataset.size = 100
h.set.speed <- c(0.2,0.5, 0.75,1.0,2.0 ,3.0, 5.0, 8.0) 

plot.set <- c("NPEM", "GP")
results.array.speed <- array(NA, dim = c(n.sample.sims, length(h.set.speed), length(plot.set), length(mu.set)))
results.array.likelihood <- array(NA, dim = c(n.sample.sims, length(h.set.speed), length(plot.set), length(mu.set)))

hyper.param.idx <- expand.grid(1:length(h.set.speed),1:length(plot.set), 1:length(mu.set))
R_bins = 1000



A.model <- A.matrix.compute(R_bins = 1000, N = 30, 
                            ker = ker.model, 
                            h = h.model, numeric.points = 400)
for(l in 1:length(h.set.speed)){
  h.tmp <- h.set.speed[l]
  A.model <- A.matrix.compute(R_bins = 1000, N = 30, 
                              ker = ker.model, 
                              h = h.tmp, numeric.points = 400)
  print(paste0("Bandwidth: ", l, " of ", length(h.set.speed)))
  for(i in 1:n.sample.sims){
    latent.true.samp <- rbeta(dataset.size, shape1 = beta.par.1, shape2 = beta.par.2)
    
    y.true.samp <- lapply(latent.true.samp, function(x){
      i = 0:N
      assumption.weights <- ker.model((i - N*x)/h.tmp)
      out <- sample(0:N, size = 1, replace = TRUE, prob = assumption.weights)
      
    })
    
    y.true.samp <- unlist(y.true.samp)
    
    y.true.frame <- matrix(data = y.true.samp, ncol = 1, byrow = TRUE)
    
    
    train.p.hat <- c()
    for(y in 0:(N)){
      prop.tmp <- mean(y.true.frame[,1] == y)
      train.p.hat <- c(train.p.hat, prop.tmp)
    }
    
    print(paste0("Dataset Sample: ", i, " of ", n.sample.sims))
    
    for(k in 1:length(mu.set)){
      mu.tmp <- mu.set[k]
      
      # time this
      # uniform initialization
      tau <- seq(0,1,length.out = R_bins)
      latent.trait.quantiles.init <- tau.init
      em.time <- system.time(em.quantiles <- nonpar.em.fit(p.hat = train.p.hat, tau = tau, 
                                    latent.trait.quantiles = latent.trait.quantiles.init, 
                                    N = 30, ker = ker.model, 
                                    J = 1000, h = h.tmp, mu = mu.tmp))

      #time this 
      gp.time <- system.time(model.estimate.binned <- numeric.latent.fit(p.hat = train.p.hat, 
                                                    A.matrix = A.model, 
                                                    mu = mu.tmp, 
                                                    show.plot = F))
      results.array.speed[i,l,1,k] <- em.time[3]
      results.array.speed[i,l,2,k] <- gp.time[3]
      
      em.p.ma <- p.ma.from.latent(tau = tau,
                                  latent.trait.quantiles = em.quantiles,
                                  ker = ker.model,h = h.tmp,N = N,
                                  mu = mu.tmp,
                                  numeric.points = 100)
      
      em.lik <- compute.likelihood(p.hat = train.p.hat, 
                                   p.ma = em.p.ma, tau = tau, 
                                   latent.trait.quantiles = em.quantiles,
                                   N = N, mu = mu.tmp)
      
      # Defines the quantiles of the latent bins 
      tau.bins <- tau.set(model.estimate.binned$latent)
      gp.p.ma <- as.vector(model.estimate.binned$observed)
      
      gp.lik <- compute.likelihood(p.hat = train.p.hat, 
                                   p.ma = gp.p.ma, tau = tau.bins, 
                                   latent.trait.quantiles = latent.trait.quantiles.init,
                                   N = N, mu = mu.tmp)
      
      results.array.likelihood[i,l,1,k] <- em.lik
      results.array.likelihood[i,l,2,k] <- gp.lik
      
    }
  }
}
  


#saveRDS(results.array.speed,paste0("data/speed_test_sims_100.rds"))
#saveRDS(results.array.likelihood,paste0("data/speed_test_likelihood_sims_100.rds"))
results.array.speed = readRDS("data/speed_test_sims_100.rds")
results.array.likelihood = readRDS("data/speed_test_likelihood_sims_100.rds")

mean.res.speed <- array(0, dim = c(length(h.set.speed), length(plot.set), length(mu.set)))

for(i in 1:n.sample.sims){
  mean.res.speed <- mean.res.speed + results.array.speed[i, , , ]
}
mean.res.speed <- mean.res.speed/n.sample.sims


sd.res.speed.tmp <- array(0, dim = c(n.sample.sims, length(h.set.speed), length(plot.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.speed.tmp[i, , , ]<-  (results.array.speed[i, , , ] - mean.res.speed)^2 
}


sd.res.speed <- array(0, dim = c(length(h.set.speed), length(plot.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.speed <- sd.res.speed + sd.res.speed.tmp[i, , , ] 
  
}

sd.res.speed <- sqrt(sd.res.speed/(n.sample.sims -1))


# standard deviation of the mean
sd.res.speed <- sd.res.speed/sqrt(n.sample.sims)



mean.res.likelihood <- array(0, dim = c(length(h.set.speed), length(plot.set), length(mu.set)))

for(i in 1:n.sample.sims){
  mean.res.likelihood <- mean.res.likelihood + results.array.likelihood[i, , , ]
}
mean.res.likelihood <- mean.res.likelihood/n.sample.sims


sd.res.likelihood.tmp <- array(0, dim = c(n.sample.sims, length(h.set.speed), length(plot.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.likelihood.tmp[i, , , ]<-  (results.array.likelihood[i, , , ] - mean.res.likelihood)^2 
}


sd.res.likelihood <- array(0, dim = c(length(h.set.speed), length(plot.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.likelihood <- sd.res.likelihood + sd.res.likelihood.tmp[i, , , ] 
  
}

sd.res.likelihood <- sqrt(sd.res.likelihood/(n.sample.sims -1))


# standard deviation of the mean
sd.res.likelihood <- sd.res.likelihood/sqrt(n.sample.sims)



### Relative Difference of likelihood: 


mean.res.diff.likelihood <- array(0, dim = c(length(h.set.speed), length(mu.set)))

for(i in 1:n.sample.sims){
  mean.res.diff.likelihood <- mean.res.diff.likelihood + (results.array.likelihood[i, ,1, ] - results.array.likelihood[i, ,2, ])/(abs(results.array.likelihood[i, ,1, ]))
}
mean.res.diff.likelihood <- mean.res.diff.likelihood/n.sample.sims


sd.res.diff.likelihood.tmp <- array(0, dim = c(n.sample.sims, length(h.set.speed), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.diff.likelihood.tmp[i, , ] <-  ((results.array.likelihood[i, ,1, ] - results.array.likelihood[i, ,2, ])/(abs(results.array.likelihood[i, ,1, ])) - mean.res.diff.likelihood)^2 
}


sd.res.diff.likelihood <- array(0, dim = c(length(h.set.speed), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.diff.likelihood <- sd.res.diff.likelihood + sd.res.diff.likelihood.tmp[i, , ] 
  
}

sd.res.diff.likelihood <- sqrt(sd.res.diff.likelihood/(n.sample.sims -1))


# standard deviation of the mean
sd.res.diff.likelihood <- sd.res.diff.likelihood/sqrt(n.sample.sims)

mus <- round(rep(mu.set, each = length(h.set.speed)), 3)
hs <- rep(h.set.speed, times = length(mu.set))

npem.speed.frame <- data.frame(h = hs, mu = mus, time = as.vector(mean.res.speed[,1,]), sd = as.vector(sd.res.speed[,1,]) )
npem.likelihood.frame <- data.frame(h = hs, mu = mus, likelihood = as.vector(mean.res.likelihood[,1,]), sd = as.vector(sd.res.likelihood[,1,]))
gp.speed.frame <- data.frame(h = hs, mu = mus, time = as.vector(mean.res.speed[,2,]), sd = as.vector(sd.res.speed[,2,]) )
gp.likelihood.frame <- data.frame(h = hs, mu = mus, likelihood = as.vector(mean.res.likelihood[,2,]), sd = as.vector(sd.res.likelihood[,2,]))
diff.likelihood.frame <- data.frame(h = hs, mu = mus, likelihood = as.vector(mean.res.diff.likelihood[,]), sd = as.vector(sd.res.diff.likelihood[,]))


npem.speed.frame$mu <- factor(npem.speed.frame$mu)
npem.likelihood.frame$mu <- factor(npem.likelihood.frame$mu)
gp.speed.frame$mu <- factor(gp.speed.frame$mu)
gp.likelihood.frame$mu <- factor(gp.likelihood.frame$mu)
diff.likelihood.frame$mu <- factor(diff.likelihood.frame$mu)



y.min.like = min(c(mean.res.likelihood[,,]))
y.max.like= max(c(mean.res.likelihood[,,]))

npem.speed.plot <- ggplot(data = npem.speed.frame, aes(x = log(h), y = time, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=time-2*sd, ymax=time+2*sd), width=.2,
                              position=position_dodge(0.05)) + ggtitle("NPEM Speed") + ylab("Fit Time (s)") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  
npem.likelihood.plot <- ggplot(data = npem.likelihood.frame, aes(x = log(h), y = likelihood, color = mu, group = mu)) + geom_line() + geom_errorbar(aes(ymin=likelihood-2*sd, ymax=likelihood+2*sd), width=.2,
                                                                                                                            position=position_dodge(0.05)) + ggtitle("NPEM Likelihood") + coord_cartesian(ylim = c(y.min.like - 0.05, y.max.like + 0.05))
gp.speed.plot <- ggplot(data = gp.speed.frame, aes(x = log(h), y = time, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=time-2*sd, ymax=time+2*sd), width=.2,
                              position=position_dodge(0.05)) + ggtitle("GP Speed") + ylab("Fit Time (s)") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

gp.likelihood.plot <- ggplot(data = gp.likelihood.frame, aes(x = log(h), y = likelihood, color = mu, group = mu)) + geom_line() + geom_errorbar(aes(ymin=likelihood-2*sd, ymax=likelihood+2*sd), width=.2,
                                                                                                                                               position=position_dodge(0.05)) + ggtitle("GP Likelihood") + coord_cartesian(ylim = c(y.min.like - 0.05, y.max.like + 0.05))
diff.likelihood.plot <- ggplot(data = diff.likelihood.frame, aes(x = log(h), y = likelihood, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=likelihood-2*sd, ymax=likelihood+2*sd), width=.2,
                              position=position_dodge(0.05)) + ggtitle("Log-Likelihood Difference (NPEM - GP)/|NPEM|") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  




grid.arrange(arrangeGrob(npem.speed.plot, gp.speed.plot, ncol = 2),                             
             diff.likelihood.plot, 
             nrow = 2)   


png(filename = "plots/SpeedTest.png",
    width = png.width, height = png.height)


grid.arrange(arrangeGrob(npem.speed.plot, gp.speed.plot, ncol = 2),                             
             diff.likelihood.plot, 
             nrow = 2)   


# Close the pdf file
dev.off() 




#########################################


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



######## SIMULATIONS ######### 
n.sample.sims = 50 
h.set<- c(0.2,0.5, 0.75,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,6.5,8.0) 

ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 3)))


results.array <- array(NA, dim = c(n.sample.sims, length(h.set), length(ker.set), length(mu.set)))
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))


# compact support kernels must have h >= 1

idx1 <- hyper.param.idx[,1]  %in% which(h.set < 1)
idx2 <- hyper.param.idx[,2]  %in% c(3,4)

hyper.param.idx <- hyper.param.idx[!(idx1 & idx2),]

for(i in 1:n.sample.sims){
  latent.true.samp <- rbeta(dataset.size, shape1 = beta.par.1, shape2 = beta.par.2)
  
  y.true.samp <- lapply(latent.true.samp, function(x){
    i = 0:N
    assumption.weights <- ker.model((i - N*x)/h.model)
    out <- sample(0:N, size = 2, replace = TRUE, prob = assumption.weights)
    
  })
  
  y.true.samp <- unlist(y.true.samp)
  
  y.true.frame <- matrix(data = y.true.samp, ncol = 2, byrow = TRUE)
  
  
  train.p.hat <- c()
  for(y in 0:(N)){
    prop.tmp <- mean(y.true.frame[,1] == y)
    train.p.hat <- c(train.p.hat, prop.tmp)
  }
  print(paste0("Dataset Sample: ", i, " of ", n.sample.sims))
  
  res.list <- lapply(1:nrow(hyper.param.idx), function(x){
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

n.sample.sims <- nrow(results.array)

mean.res <- array(0, dim = c(length(h.set), length(ker.set), length(mu.set)))

for(i in 1:n.sample.sims){
  mean.res <- mean.res + results.array[i, , , ]
}
mean.res <- mean.res/n.sample.sims


sd.res.tmp <- array(0, dim = c(n.sample.sims, length(h.set), length(ker.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.tmp[i, , , ]<-  (results.array[i, , , ] - mean.res)^2 
}


sd.res <- array(0, dim = c(length(h.set), length(ker.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res <- sd.res + sd.res.tmp[i, , , ] 
  
}

sd.res <- sqrt(sd.res/(n.sample.sims -1))


# standard deviation of the mean
sd.res <- sd.res/sqrt(n.sample.sims)

mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(mu.set))

gauss.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(mean.res[,1,]), total.variation.sd = as.vector(sd.res[,1,]) )
exp.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(mean.res[,2,]), total.variation.sd = as.vector(sd.res[,2,]))
triangle.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(mean.res[,3,]), total.variation.sd = as.vector(sd.res[,3,]))
epanechnikov.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(mean.res[,4,]), total.variation.sd = as.vector(sd.res[,4,]))


# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Gaussian Kernel") + 
  ylab("Total Variation") + 
  geom_vline(xintercept = h.model) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Laplace Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Triangle Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))

epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2,position=position_dodge(0.05)) + ggtitle("Epanechnikov Kernel") + 
  ylab("Total Variation") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  + 
  coord_cartesian(ylim = c(.1, .7))


gauss.plot
exp.plot
triangle.plot
epanechnikov.plot




png(filename = "plots/intrinsic_variability_sim1_large.png",
    width = png.width, height = png.height)

ggarrange(gauss.plot, exp.plot, 
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 

png(filename = "plots/intrinsic_variability_sim1.png",
    width = 2*png.width, height = png.height)

ggarrange(gauss.plot, exp.plot, epanechnikov.plot,
          ncol=3, nrow=1, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 

#saveRDS(results.array,paste0("data/intrinsic_variability_sims_100.rds"))





######### Variation as a Function of mu within the proper model 
set.seed(1)
ker.correct.model <- ker.model
h.correct.model <- 2  #h.model
A.correct.model <- A.matrix.compute(R_bins = 1000, N = 30, 
                                    ker = ker.correct.model, 
                                    h = h.correct.model, 
                                    numeric.points = 400)

A.two.sample <- A.two.obs.tensor.compute(R_bins = 1000, N = 30, 
                                         ker = ker.correct.model, 
                                         h = h.correct.model, 
                                         numeric.points = 400)


beta.par.1 <- 12 
beta.par.2 <- 5


dataset.size = 100
n.sample.sims <- 200
long.mu.set <- c(0,exp(seq(log(0.001),log(.1), length.out = 15)))

results.array.long.mu <- array(NA, dim = c(n.sample.sims, length(long.mu.set)))
hyper.param.idx <- expand.grid(1:length(long.mu.set))


for(i in 1:n.sample.sims){
  latent.true.samp <- rbeta(dataset.size, shape1 = beta.par.1, shape2 = beta.par.2)
  
  y.true.samp <- lapply(latent.true.samp, function(x){
    i = 0:N
    assumption.weights <- ker.correct.model((i - N*x)/h.correct.model)
    out <- sample(0:N, size = 2, replace = TRUE, prob = assumption.weights)
    
  })
  
  y.true.samp <- unlist(y.true.samp)
  
  y.true.frame <- matrix(data = y.true.samp, ncol = 2, byrow = TRUE)
  
  
  train.p.hat <- c()
  for(y in 0:(N)){
    prop.tmp <- mean(y.true.frame[,1] == y)
    train.p.hat <- c(train.p.hat, prop.tmp)
  }
  
  val.p.hat <- c()
  for(y in 0:(N)){
    prop.tmp <- mean(y.true.frame[,2] == y)
    val.p.hat <- c(val.p.hat, prop.tmp)
  }
  
  
  print(paste0("Dataset Sample: ", i, " of ", n.sample.sims))
  
  res.list <- lapply(1:nrow(hyper.param.idx), function(x){
    l = hyper.param.idx[x,1]
    #h.tmp <- h.correct.model
    #ker.tmp <- ker.correct.model
    mu.tmp <- long.mu.set[l]
    
    model.estimate <- numeric.latent.fit(p.hat = train.p.hat, 
                                         A.matrix = A.correct.model, 
                                         mu = mu.tmp, show.plot = F)
    
    
    latent.mix.list <- list()
    
    for (m in 1:nrow(y.true.frame)) {
      latent.mix.list[[m]] <- model.estimate$latent
    }
    
    
    log.like <- two.samp.log.likelihood(y.true.frame, A.two.sample, latent.mix.list)
    
    
    res.tmp <- -log.like
    return(res.tmp)
  })
  res.set <- unlist(res.list)
  for(q in 1:nrow(hyper.param.idx)){
    results.array.long.mu[i,q] <- res.set[q]
  }
  
}


#saveRDS(results.array.long.mu,paste0("data/correct_model_regularization_sims.rds"))
results.array.long.mu <- readRDS("data/correct_model_regularization_sims.rds")

n.sample.sims <- nrow(results.array.long.mu)

mean.res.long.mu <- array(0, dim = c(length(long.mu.set)))

for(i in 1:n.sample.sims){
  mean.res.long.mu <- mean.res.long.mu + as.numeric(results.array.long.mu[i, ])
}

mean.res.long.mu <- mean.res.long.mu/n.sample.sims


sd.res.tmp <- array(0, dim = c(n.sample.sims, length(long.mu.set)))
for(i in 1:n.sample.sims){
  sd.res.tmp[i,] <-  as.numeric((results.array.long.mu[i,] - mean.res.long.mu)^2)
}


sd.res.long.mu <- array(0, dim = c(length(long.mu.set)))
for(i in 1:n.sample.sims){
  sd.res.long.mu <- sd.res.long.mu + sd.res.tmp[i, ] 
  
}

sd.res.long.mu <- sqrt(sd.res.long.mu/(n.sample.sims -1))
sd.of.mean <-  sd.res.long.mu/ sqrt(n.sample.sims)



long.mu.frame <- data.frame(log.mu = log(long.mu.set + 0.001), total.variation = mean.res.long.mu, total.variation.sd = sd.of.mean)
long.mu.plot <- ggplot(data = long.mu.frame, aes(x = log.mu, y = total.variation)) + geom_line() + 
  geom_errorbar(aes(ymin=total.variation-2*total.variation.sd, ymax=total.variation+2*total.variation.sd), 
                width=.05, position=position_dodge(0.05)) + ylab("Negative Two Sample Log-Likelihood") + 
  ggtitle("Latent: Beta(12,5) and MKM(Gaussian,h = 2)") + xlab("log(\u00b5 + 0.001)") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  
long.mu.plot



print(paste0("Optimal Mu Value: ", long.mu.set[which.min(mean.res.long.mu)]))
# 0.0193





png(filename = "plots/mu_selection_sim1.png",
    width = png.width, height = png.height)

long.mu.plot

# Close the pdf file
dev.off() 





### repeat with a second models 
######### Variation as a Function of mu within the proper model 
set.seed(1)
ker.correct.model <- exponential_kernel
h.correct.model <- 1  #h.model
A.correct.model <- A.matrix.compute(R_bins = 1000, N = 30, 
                                    ker = ker.correct.model, 
                                    h = h.correct.model, 
                                    numeric.points = 400)

A.two.sample <- A.two.obs.tensor.compute(R_bins = 1000, N = 30, 
                                         ker = ker.correct.model, 
                                         h = h.correct.model, 
                                         numeric.points = 400)


beta.par.1 <- 6 
beta.par.2 <- 6



dataset.size = 100
n.sample.sims <- 500
long.mu.set <- c(0,exp(seq(log(0.001),log(.1), length.out = 15)))

results.array.long.mu <- array(NA, dim = c(n.sample.sims, length(long.mu.set)))
hyper.param.idx <- expand.grid(1:length(long.mu.set))


for(i in 1:n.sample.sims){
  latent.true.samp <- rbeta(dataset.size, shape1 = beta.par.1, shape2 = beta.par.2)
  
  y.true.samp <- lapply(latent.true.samp, function(x){
    i = 0:N
    assumption.weights <- ker.correct.model((i - N*x)/h.correct.model)
    out <- sample(0:N, size = 2, replace = TRUE, prob = assumption.weights)
    
  })
  
  y.true.samp <- unlist(y.true.samp)
  
  y.true.frame <- matrix(data = y.true.samp, ncol = 2, byrow = TRUE)
  
  
  train.p.hat <- c()
  for(y in 0:(N)){
    prop.tmp <- mean(y.true.frame[,1] == y)
    train.p.hat <- c(train.p.hat, prop.tmp)
  }
  
  val.p.hat <- c()
  for(y in 0:(N)){
    prop.tmp <- mean(y.true.frame[,2] == y)
    val.p.hat <- c(val.p.hat, prop.tmp)
  }
  
  
  print(paste0("Dataset Sample: ", i, " of ", n.sample.sims))
  
  res.list <- lapply(1:nrow(hyper.param.idx), function(x){
    l = hyper.param.idx[x,1]
    #h.tmp <- h.correct.model
    #ker.tmp <- ker.correct.model
    mu.tmp <- long.mu.set[l]
    
    model.estimate <- numeric.latent.fit(p.hat = train.p.hat, 
                                         A.matrix = A.correct.model, 
                                         mu = mu.tmp, show.plot = F)
    
    
    latent.mix.list <- list()
    
    for (m in 1:nrow(y.true.frame)) {
      latent.mix.list[[m]] <- model.estimate$latent
    }
    
    
    log.like <- two.samp.log.likelihood(y.true.frame, A.two.sample, latent.mix.list)
    
    
    res.tmp <- -log.like
    return(res.tmp)
  })
  res.set <- unlist(res.list)
  for(q in 1:nrow(hyper.param.idx)){
    results.array.long.mu[i,q] <- res.set[q]
  }
  
}

mean.res.long.mu <- array(0, dim = c(length(long.mu.set)))

for(i in 1:n.sample.sims){
  mean.res.long.mu <- mean.res.long.mu + results.array.long.mu[i, ]
}

mean.res.long.mu <- mean.res.long.mu/n.sample.sims


sd.res.tmp <- array(0, dim = c(n.sample.sims, length(long.mu.set)))
for(i in 1:n.sample.sims){
  sd.res.tmp[i,]<-  (results.array.long.mu[i,] - mean.res.long.mu)^2 
}


sd.res.long.mu <- array(0, dim = c(length(long.mu.set)))
for(i in 1:n.sample.sims){
  sd.res.long.mu <- sd.res.long.mu + sd.res.tmp[i, ] 
  
}

sd.res.long.mu <- sqrt(sd.res.long.mu/(n.sample.sims -1))


saveRDS(results.array.long.mu,paste0("data/correct_model_regularization_sims_2.rds"))


sd.of.mean <-  sd.res.long.mu/ sqrt(n.sample.sims)

long.mu.frame <- data.frame(log.mu = log(long.mu.set + 0.001), total.variation = mean.res.long.mu, total.variation.sd = sd.of.mean)
long.mu.plot <- ggplot(data = long.mu.frame, aes(x = log.mu, y = total.variation)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-2*total.variation.sd, 
                                  ymax=total.variation+2*total.variation.sd), 
                              width=.05,position=position_dodge(0.05)) + 
  ylab("Negative Two Sample Log-Likelihood")+ ggtitle("Latent: Beta(6,6) and MKM(Laplace,h = 1)") + xlab("log(\u00b5 + 0.001)") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

long.mu.plot




print(paste0("Optimal Mu Value: ", long.mu.set[which.min(mean.res.long.mu)]))
# 0.01390 
















############################################################
############################################################

# two scores generated under the following model 
beta.par.1y <- 12 
beta.par.2y <- 5
h.model.y <- 2
ker.model.y <- gaussian_kernel

beta.par.1z <- 6 
beta.par.2z <- 6
h.model.z <- 1
ker.model.z <- exponential_kernel

n.train.val <- 100
#n.test <- 500 # for monte carlo approximation of test set error 

mu.y <- 0.01930
mu.z <- 0.01390 
# Assume we have already computed an optimal h ker and mu from the first model 
Ny = 30
Nz = 30

num.bins <- 300


A.model.y <- A.matrix.compute(R_bins = num.bins, N = Ny, 
                              ker = ker.model.y, h = h.model.y, numeric.points = 400)


h.set<- c(0.05,0.5,1.0,1.5,2.0,2.5,4.0) 
#h.set<- c(0.05,0.5,1.0,1.5,2.0) 
ker.set <- list(gaussian_kernel, exponential_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(mu.z), length.out = 3)), exp(seq(log(mu.y),log(.3), length.out = 3)))

#mu.set <- c(0,0.1)







# Pre computing A matrices
A.matrix.set <- list()
for(j in 1:length(h.set)){
  h.tmp <- h.set[j]
  A.matrix.set.tmp <- list()
  for(k in 1:length(ker.set)){
    
    ker.tmp <- ker.set[[k]]
    
    A.mat <- A.matrix.compute(R_bins = num.bins, N = 30, ker = ker.tmp, h = h.tmp, numeric.points = 400)
    A.matrix.set.tmp[[k]] <- A.mat 
    
  }
  A.matrix.set[[j]] <- A.matrix.set.tmp
}



######## SIMULATIONS ######### 
n.sample.sims = 10 
results.array <- array(NA, dim = c(n.sample.sims, length(h.set), length(ker.set), length(mu.set)))
results.array.ml <- array(NA, dim = c(n.sample.sims, length(h.set), length(ker.set)))
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))
hyper.param.idx.ml <- expand.grid(1:length(h.set),1:length(ker.set))
  

# compact support kernels must have h >= 1

idx1 <- hyper.param.idx[,1]  %in% which(h.set < 1)
idx2 <- hyper.param.idx[,2]  %in% c(3,4)

hyper.param.idx <- hyper.param.idx[!(idx1 & idx2),]


### population test set 


grid.size <- 10000
latent.grid <- seq(0,1, length.out = grid.size)

test.joint.dist <- lapply(latent.grid, function(x){
  latent.y <- qbeta(x, shape1 = beta.par.1y, shape2 = beta.par.2y)
  latent.z <- qbeta(x, shape1 = beta.par.1z, shape2 = beta.par.2z)
  iy = 0:Ny
  iz = 0:Nz
  assumption.weights.y <- ker.model.y((iy - Ny*latent.y)/h.model.y)
  assumption.weights.z <- ker.model.z((iz - Nz*latent.z)/h.model.z)
  
  prob.y <- assumption.weights.y/sum(assumption.weights.y)
  prob.z <- assumption.weights.z/sum(assumption.weights.z)
  
  p.yz <- prob.y %*% t(prob.z)
  out <- p.yz
  
})

p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
for(i in 1:length(test.joint.dist)){
  p.yz <- p.yz + test.joint.dist[[i]]
}
p.yz <- p.yz/grid.size
test.grid <- expand.grid(0:Ny, 0:Nz)


library(plot.matrix)
plot(p.yz)



for(i in 1:n.sample.sims){
  
  latent.uniform <- runif(n.train.val)

  pair.true.samp <- lapply(latent.uniform, function(x){
    latent.y <- qbeta(x, shape1 = beta.par.1y, shape2 = beta.par.2y)
    latent.z <- qbeta(x, shape1 = beta.par.1z, shape2 = beta.par.2z)
    iy = 0:Ny
    iz = 0:Nz
    assumption.weights.y <- ker.model.y((iy - Ny*latent.y)/h.model.y)
    assumption.weights.z <- ker.model.z((iz - Nz*latent.z)/h.model.z)
    
    out.y <- sample(0:Ny, size = 1, replace = TRUE, prob = assumption.weights.y)
    out.z <- sample(0:Nz, size = 1, replace = TRUE, prob = assumption.weights.z)
    out <- c(out.y, out.z)
    
  })
  
  pair.true.samp <- unlist(pair.true.samp)
  pair.true.samp <- matrix(data = pair.true.samp, ncol = 2, byrow = TRUE)
  
  
  train.p.hat.y <- c()
  for(y in 0:(Ny)){
    prop.tmp <- mean(pair.true.samp[,1] == y)
    train.p.hat.y <- c(train.p.hat.y, prop.tmp)
  }
  
  ### Fixed Model Estimate of y under correctly learned model 
  model.estimate.y <- numeric.latent.fit(p.hat = train.p.hat.y, 
                                       A.matrix = A.model.y, 
                                       mu = mu.y, show.plot = T)
  
  model.estimate.y.ml <- numeric.latent.fit(p.hat = train.p.hat.y, 
                                           A.matrix = A.model.y, 
                                           mu = 0, show.plot = T)
  
  train.p.hat.z <- c()
  for(z in 0:(Nz)){
    prop.tmp <- mean(pair.true.samp[,2] == z)
    train.p.hat.z <- c(train.p.hat.z, prop.tmp)
  }
  
  print(paste0("Dataset Sample: ", i, " of ", n.sample.sims))
  
  res <- sapply(1:nrow(hyper.param.idx), function(x){
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    l = hyper.param.idx[x,3]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- mu.set[l]
    A.matrix <- A.matrix.set[[j]][[k]]
    model.estimate.z <- numeric.latent.fit(p.hat = train.p.hat.z, 
                                           A.matrix = A.matrix, 
                                           mu = mu.tmp, show.plot = F)
    
    latent.mix.list.y <- list()
    model.observed.list.y <- list()
    latent.mix.list.z <- list()
    model.observed.list.z <- list()
    
    
    for (m in 1:nrow(test.grid)) {
      latent.mix.list.y[[m]] <- model.estimate.y$latent
      model.observed.list.y[[m]] <- model.estimate.y$observed
      latent.mix.list.z[[m]] <- model.estimate.z$latent
      model.observed.list.z[[m]] <- model.estimate.z$observed
      
    }
    
    res.tmp <- score.conversion(test.pairs = test.grid, latent.mix.list.y = latent.mix.list.y,
                                model.observed.list.y = model.observed.list.y,
                                latent.mix.list.z = latent.mix.list.z, 
                                model.observed.list.z = model.observed.list.z, 
                                Ny = Ny, ker.y = ker.model.y, hy = h.model.y,
                                Nz = Nz, ker.z = ker.tmp, hz = h.tmp, joint.prob = p.yz,
                                n.samp = 4, method = "CrossEntropy")
    
    
    cat(paste0("Param Set: ", x, "/", nrow(hyper.param.idx)), end = "\r")
    return(res.tmp)
  })
  cat( end = "\n")
  res.set <- res
  
  res.set <- res
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    l = hyper.param.idx[q,3]
    results.array[i,j,k,l] <- res.set[q]
  }
  
  
  cat(paste0("ML conversion"), end = "\n")
  
  ##### ML conversion 
  
  res <- sapply(1:nrow(hyper.param.idx.ml), function(x){
    j = hyper.param.idx.ml[x,1]
    k = hyper.param.idx.ml[x,2]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    
    A.matrix <- A.matrix.set[[j]][[k]]
    model.estimate.z.ml <- numeric.latent.fit(p.hat = train.p.hat.z, 
                                           A.matrix = A.matrix, 
                                           mu = 0, show.plot = F)
    
    latent.mix.list.y.ml <- list()
    model.observed.list.y.ml <- list()
    latent.mix.list.z.ml <- list()
    model.observed.list.z.ml <- list()
    
    # 
    for (m in 1:nrow(test.grid)) {
      latent.mix.list.y.ml[[m]] <- model.estimate.y.ml$latent
      model.observed.list.y.ml[[m]] <- model.estimate.y.ml$observed
      latent.mix.list.z.ml[[m]] <- model.estimate.z.ml$latent
      model.observed.list.z.ml[[m]] <- model.estimate.z.ml$observed
      
    }
    
    res.tmp <- score.conversion(test.pairs = test.grid, 
                                latent.mix.list.y = latent.mix.list.y.ml,
                                model.observed.list.y = model.observed.list.y.ml,
                                latent.mix.list.z = latent.mix.list.z.ml, 
                                model.observed.list.z = model.observed.list.z.ml, 
                                Ny = Ny, ker.y = ker.model.y, hy = h.model.y,
                                Nz = Nz, ker.z = ker.tmp, hz = h.tmp, joint.prob = p.yz,
                                n.samp = 4, method = "CrossEntropy")
    
    
    cat(paste0("Param Set: ", x, "/", nrow(hyper.param.idx.ml)), end = "\r")
    return(res.tmp)
  })
  cat( end = "\n")
  res.set <- res
  for(q in 1:nrow(hyper.param.idx.ml)){
    j = hyper.param.idx.ml[q,1]
    k = hyper.param.idx.ml[q,2]
    results.array.ml[i,j,k] <- res.set[q]
  }
}


#saveRDS(results.array,paste0("data/conversion_sims_100.rds"))
#saveRDS(results.array.ml,paste0("data/conversion_ml_sims_100.rds"))

results.array <- readRDS("data/conversion_sims_100.rds")
results.array.ml <- readRDS("data/conversion_ml_sims_100.rds")

mean.res <- array(0, dim = c(length(h.set), length(ker.set), length(mu.set)))

n.sample.sims <- nrow(results.array)
for(i in 1:n.sample.sims){
  mean.res <- mean.res + results.array[i, , , ]
}
mean.res <- mean.res/n.sample.sims


sd.res.tmp <- array(0, dim = c(n.sample.sims, length(h.set), length(ker.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res.tmp[i, , , ]<-  (results.array[i, , , ] - mean.res)^2 
}


sd.res <- array(0, dim = c(length(h.set), length(ker.set), length(mu.set)))
for(i in 1:n.sample.sims){
  sd.res <- sd.res + sd.res.tmp[i, , , ] 
  
}

sd.res <- sqrt(sd.res/(n.sample.sims -1))




mean.res.ml <- array(0, dim = c(length(h.set), length(ker.set)))

for(i in 1:n.sample.sims){
  mean.res.ml <- mean.res.ml + results.array.ml[i, , ]
}
mean.res.ml <- mean.res.ml/n.sample.sims


sd.res.tmp.ml <- array(0, dim = c(n.sample.sims, length(h.set), length(ker.set)))
for(i in 1:n.sample.sims){
  sd.res.tmp.ml[i, , ]<-  (results.array.ml[i, , ] - mean.res.ml)^2 
}


sd.res.ml <- array(0, dim = c(length(h.set), length(ker.set)))
for(i in 1:n.sample.sims){
  sd.res.ml <- sd.res.ml + sd.res.tmp.ml[i, , ] 
  
}

sd.res.ml <- sqrt(sd.res.ml/(n.sample.sims -1))




library(Polychrome)

nb.cols <- 4
P4 = createPalette(nb.cols,  c("#ff0000", "#00ff00", "#0000ff"))

colfunc <- colorRampPalette(c("black", "white"))
colfunc(10)



palette(P4)


library(RColorBrewer)
display.brewer.all()

P5 <- brewer.pal(6, "Blues")
P5 <- P5[2:6]
P6 <- c(P5,"red")
mus <- round(rep(mu.set, each = length(h.set)), 4)
mus <- as.character(mus)
mus <- c(mus, rep("unregularized", length(h.set)))
hs <- rep(h.set, times = length(mu.set) + 1)

gauss.frame <- data.frame(h = hs, mu = mus, total.variation = c(as.vector(mean.res[,1,]),as.vector(mean.res.ml[,1])), total.variation.sd = c(as.vector(sd.res[,1,]),as.vector(sd.res.ml[,1])) )
exp.frame <- data.frame(h = hs, mu = mus, total.variation = c(as.vector(mean.res[,2,]),as.vector(mean.res.ml[,2])), total.variation.sd = c(as.vector(sd.res[,2,]),as.vector(sd.res.ml[,2])))
# triangle.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(mean.res[,3,]), total.variation.sd = as.vector(sd.res[,3,]))
# epanechnikov.frame <- data.frame(h = hs, mu = mus, total.variation = as.vector(mean.res[,4,]), total.variation.sd = as.vector(sd.res[,4,]))


#Narrowing Colour Choices 
gauss.frame <- gauss.frame[gauss.frame$mu %in% c("0","0.0037", "0.0193", "0.0761","0.3","unregularized"), ]
exp.frame <- exp.frame[exp.frame$mu %in% c("0","0.0037", "0.0193", "0.0761","0.3","unregularized"), ]


gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
# triangle.frame$mu <- factor(triangle.frame$mu)
# epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)



# model selection index, exp frame 31 
gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2, position=position_dodge(0.0)) + ggtitle("Gaussian Kernel") + 
  coord_cartesian(xlim = c(0,4.1), ylim = c(2.5, 3)) + ylab("Population Cross Entropy")  

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), 
                              width=.2, position=position_dodge(0.0)) + ggtitle("Laplace Kernel") + 
  coord_cartesian(xlim = c(0,4.1), ylim = c(2.5, 3)) + ylab("Population Cross Entropy")  + 
  annotate("point", x = exp.frame$h[17], y = exp.frame$total.variation[17], colour = "blue", shape = "x", size = 8) #+ geom_point(size = 0.6) +  geom_point(aes(x = exp.frame$h[31], y = exp.frame$total.variation[31], shape = "+", color = "black"))
# triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), width=.2,
#                                                                                                                                       position=position_dodge(0.05)) + ggtitle("Triangle Kernel")
# epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + geom_errorbar(aes(ymin=total.variation-3*total.variation.sd, ymax=total.variation+3*total.variation.sd), width=.2,
#                                                                                                                                               position=position_dodge(0.05)) + ggtitle("Epanechnikov Kernel")

# argmin is the gaussian kernel with mu = 0.01




gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line()  + ggtitle("Gaussian Kernel")  + coord_cartesian(xlim = c(0,4.1), ylim = c(2.45, 3)) + 
  ylab("Population Cross Entropy") + scale_color_manual(values= as.vector(P6)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + 
  geom_line()  + ggtitle("Laplace Kernel") + coord_cartesian(xlim = c(0,4.1), ylim = c(2.45, 3)) + 
  ylab("Population Cross Entropy")  + 
  annotate("point", x = exp.frame$h[17], y = exp.frame$total.variation[17], colour = "blue", shape = "x", size = 8) + 
  scale_color_manual(values= as.vector(P6)) + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

# triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Triangle Kernel")
# epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Epanechnikov Kernel")






png(filename = "plots/Conversion_Simulation.png",
    width = 1.5*png.width, height = png.height)



ggarrange(gauss.plot, exp.plot,
          ncol=2, nrow=1, common.legend = TRUE, legend="right") 
# Close the pdf file
dev.off() 




############### AAIC Table 



######## SIMULATIONS ######### 

# We use 4 methods, Z-Score, Quantile Matching, NPMLE, RNPMLE 
n.sample.sims = 200 
results.array <- array(NA, dim = c(n.sample.sims, 4))


grid.size <- 10000
latent.grid <- seq(0,1, length.out = grid.size)

num.bins <- 1000
A.model.y <- A.matrix.compute(R_bins = num.bins, N = Ny, 
                              ker = ker.model.y, h = h.model.y, numeric.points = 400)

A.model.z <- A.matrix.compute(R_bins = num.bins, N = Ny, 
                              ker = ker.model.z, h = h.model.z, numeric.points = 400)



test.joint.dist <- lapply(latent.grid, function(x){
  latent.y <- qbeta(x, shape1 = beta.par.1y, shape2 = beta.par.2y)
  latent.z <- qbeta(x, shape1 = beta.par.1z, shape2 = beta.par.2z)
  iy = 0:Ny
  iz = 0:Nz
  assumption.weights.y <- ker.model.y((iy - Ny*latent.y)/h.model.y)
  assumption.weights.z <- ker.model.z((iz - Nz*latent.z)/h.model.z)
  
  prob.y <- assumption.weights.y/sum(assumption.weights.y)
  prob.z <- assumption.weights.z/sum(assumption.weights.z)
  
  p.yz <- prob.y %*% t(prob.z)
  out <- p.yz
  
})

p.yz <- matrix(data = 0, nrow = Ny + 1, ncol = Nz + 1)
for(i in 1:length(test.joint.dist)){
  p.yz <- p.yz + test.joint.dist[[i]]
}
p.yz <- p.yz/grid.size
test.grid <- expand.grid(0:Ny, 0:Nz)


library(plot.matrix)
plot(p.yz)


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

quantile.conversion.prob <- function(y,quantile_map_yz, decay_rate, Nz){
  z.pred <- quantile_map_yz[y + 1]
  
  z.scale <- 0:Nz
  probs <- exp(-decay_rate*abs(z.pred - z.scale))
  
  out <- probs/sum(probs)
  return(out)
}

for(i in 1:n.sample.sims){
  
  latent.uniform <- runif(n.train.val)
  
  pair.true.samp <- lapply(latent.uniform, function(x){
    latent.y <- qbeta(x, shape1 = beta.par.1y, shape2 = beta.par.2y)
    latent.z <- qbeta(x, shape1 = beta.par.1z, shape2 = beta.par.2z)
    iy = 0:Ny
    iz = 0:Nz
    assumption.weights.y <- ker.model.y((iy - Ny*latent.y)/h.model.y)
    assumption.weights.z <- ker.model.z((iz - Nz*latent.z)/h.model.z)
    
    out.y <- sample(0:Ny, size = 1, replace = TRUE, prob = assumption.weights.y)
    out.z <- sample(0:Nz, size = 1, replace = TRUE, prob = assumption.weights.z)
    out <- c(out.y, out.z)
    
  })
  
  pair.true.samp <- unlist(pair.true.samp)
  pair.true.samp <- matrix(data = pair.true.samp, ncol = 2, byrow = TRUE)
  

  train.p.hat.y <- c()
  for(y in 0:(Ny)){
    prop.tmp <- mean(pair.true.samp[,1] == y)
    train.p.hat.y <- c(train.p.hat.y, prop.tmp)
  }

  ### Fixed Model Estimate of y under correctly learned model
  model.estimate.y <- numeric.latent.fit(p.hat = train.p.hat.y,
                                         A.matrix = A.model.y,
                                         mu = mu.y, show.plot = T)

  model.estimate.y.ml <- numeric.latent.fit(p.hat = train.p.hat.y,
                                            A.matrix = A.model.y,
                                            mu = 0, show.plot = T)

  train.p.hat.z <- c()
  for(z in 0:(Nz)){
    prop.tmp <- mean(pair.true.samp[,2] == z)
    train.p.hat.z <- c(train.p.hat.z, prop.tmp)
  }




  model.estimate.z <- numeric.latent.fit(p.hat = train.p.hat.z,
                                         A.matrix = A.model.z,
                                         mu = mu.z, show.plot = F)

  model.estimate.z.ml <- numeric.latent.fit(p.hat = train.p.hat.z,
                                         A.matrix = A.model.z,
                                         mu = 0, show.plot = F)

  latent.mix.list.y <- list()
  model.observed.list.y <- list()
  latent.mix.list.z <- list()
  model.observed.list.z <- list()


  for (m in 1:nrow(test.grid)) {
    latent.mix.list.y[[m]] <- model.estimate.y$latent
    model.observed.list.y[[m]] <- model.estimate.y$observed
    latent.mix.list.z[[m]] <- model.estimate.z$latent
    model.observed.list.z[[m]] <- model.estimate.z$observed
  }

  res <- score.conversion(test.pairs = test.grid, latent.mix.list.y = latent.mix.list.y,
                              model.observed.list.y = model.observed.list.y,
                              latent.mix.list.z = latent.mix.list.z,
                              model.observed.list.z = model.observed.list.z,
                              Ny = Ny, ker.y = ker.model.y, hy = h.model.y,
                              Nz = Nz, ker.z = ker.model.z, hz = h.model.z, joint.prob = p.yz,
                              n.samp = 4, method = "CrossEntropy")



  latent.mix.list.y.ml <- list()
  model.observed.list.y.ml <- list()
  latent.mix.list.z.ml <- list()
  model.observed.list.z.ml <- list()


  for(m in 1:nrow(test.grid)) {
    latent.mix.list.y.ml[[m]] <- model.estimate.y.ml$latent
    model.observed.list.y.ml[[m]] <- model.estimate.y.ml$observed
    latent.mix.list.z.ml[[m]] <- model.estimate.z.ml$latent
    model.observed.list.z.ml[[m]] <- model.estimate.z.ml$observed
  }


  res.ml <- score.conversion(test.pairs = test.grid, latent.mix.list.y = latent.mix.list.y.ml,
                              model.observed.list.y = model.observed.list.y.ml,
                              latent.mix.list.z = latent.mix.list.z.ml,
                              model.observed.list.z = model.observed.list.z.ml,
                              Ny = Ny, ker.y = ker.model.y, hy = h.model.y,
                              Nz = Nz, ker.z = ker.model.z, hz = h.model.z, joint.prob = p.yz,
                              n.samp = 4, method = "CrossEntropy")




  # naive method:

  muy <- mean(pair.true.samp[,1])
  sdy <- sd(pair.true.samp[,1])
  muz <- mean(pair.true.samp[,2])
  sdz <- sd(pair.true.samp[,2])



  ce.zscore = 0

  for(k in 1:nrow(test.grid)){
    score1 <- test.grid[k,1]
    score2 <- test.grid[k,2]
    z.prob = naive.conversion.prob(y = score1, muy = muy, sdy = sdy,
                                   muz = muz, sdz = sdz, Nz = Nz, n.samp = 10000)
    ce.zscore = ce.zscore -p.yz[score1 + 1,score2 + 1]*log(z.prob[score2 + 1])
    cat(paste0("Test sample ", k, "/",nrow(test.grid)), end = "\r")
  }
  cat(end = "\n")
  
  
  
  
  
  
  # quantile matching method: 
  
  # naive method: 
  
  quantile_map_yz <- c()
  for(y in 0:Ny){
    cdf.y <- mean(  pair.true.samp[,1] <= y )
    z.pred <- round(quantile(pair.true.samp[,2], cdf.y))
    quantile_map_yz[y+1] <- z.pred
  }
  
  
  
  

  
  
  ce.quant.match = 0
  
  for(k in 1:nrow(test.grid)){
    score1 <- test.grid[k,1]
    score2 <- test.grid[k,2]
    z.prob = quantile.conversion.prob(y = score1,quantile_map_yz = quantile_map_yz, decay_rate = 0.1, Nz = Nz, n.samp = 100000)
    ce.quant.match = ce.quant.match -p.yz[score1 + 1,score2 + 1]*log(z.prob[score2 + 1])
    cat(paste0("Test sample ", k, "/",nrow(test.grid)), end = "\r")
  }
  
  results.array[i,1] <- ce.zscore
  results.array[i,2] <- ce.quant.match
  results.array[i,3] <- res.ml
  results.array[i,4] <- res
    
  cat(end = "\n")
  print(paste0("Dataset Sample: ", i, " of ", n.sample.sims))
}

saveRDS(results.array,paste0("data/AAIC_sim_table.rds"))


aaic_table <- data.frame(matrix(data = NA, nrow = 4, ncol = 3))
colnames(aaic_table) <- c("Method", "Mean Cross Entropy", "SE")
aaic_table[,1] <- c("Z Score", "Quantile Matching", "Nonparametric",
                    "Regularized Nonparametric")

aaic_table[,2] <- c(mean(results.array[,1]), 
                    mean(results.array[,2]),
                    mean(results.array[,3]),
                    mean(results.array[,4]))

aaic_table[,3] <- c(sd(results.array[,1])/sqrt(n.sample.sims), 
                    sd(results.array[,2])/sqrt(n.sample.sims),
                    sd(results.array[,3])/sqrt(n.sample.sims),
                    sd(results.array[,4])/sqrt(n.sample.sims))

aaic_table[,c(2,3)] <- round(aaic_table[,c(2,3)], 3)

aaic_table
stargazer::stargazer(aaic_table, summary = NULL)

############################################################
############################################################




####### First and second order feasibility tests 
n.sample.sims = 300
dataset.size = 100
h.set <- c(0.25,0.5,1,2,4,7,10,14,20,25,30)
ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
ker.model <- gaussian_kernel
h.model <- 2
num.bins <- 300
results.array.order.1 <- array(NA, dim = c(n.sample.sims, length(h.set), length(ker.set)))
results.array.order.2 <- array(NA, dim = c(n.sample.sims, length(h.set), length(ker.set)))

hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set))

# compact support kernels must have h >= 1

idx1 <- hyper.param.idx[,1]  %in% which(h.set < 1)
idx2 <- hyper.param.idx[,2]  %in% c(3,4)

hyper.param.idx <- hyper.param.idx[!(idx1 & idx2),]


# Pre-computing A matrices 

# Pre computing A matrices
A.matrix.set.order.1 <- list()
for(j in 1:length(h.set)){
  h.tmp <- h.set[j]
  A.matrix.set.tmp <- list()
  for(k in 1:length(ker.set)){
    
    ker.tmp <- ker.set[[k]]
    
    A.mat <- A.matrix.compute(R_bins = num.bins, N = 30, ker = ker.tmp, h = h.tmp, numeric.points = 400)
    A.matrix.set.tmp[[k]] <- A.mat 
    
  }
  A.matrix.set.order.1[[j]] <- A.matrix.set.tmp
}



two.obs.mapping <- expand.grid(0:N,0:N)
A.matrix.set.order.2 <- list()
for(j in 1:length(h.set)){
  h.tmp <- h.set[j]
  A.matrix.set.tmp <- list()
  for(k in 1:length(ker.set)){
    
    ker.tmp <- ker.set[[k]]
    
    A.mat <- A.matrix.compute.two.obs(R_bins = num.bins, N = 30, ker = ker.tmp, h = h.tmp, 
                              two.obs.mapping = two.obs.mapping, numeric.points = 400)
    A.matrix.set.tmp[[k]] <- A.mat 
    
  }
  A.matrix.set.order.2[[j]] <- A.matrix.set.tmp
}



set.seed(1)
for(i in 1:n.sample.sims){
  
  latent.true.samp <- rbeta(dataset.size, shape1 = beta.par.1, shape2 = beta.par.2)
  
  y.true.samp <- lapply(latent.true.samp, function(x){
    i = 0:N
    assumption.weights <- ker.model((i - N*x)/h.model)
    out <- sample(0:N, size = 2, replace = TRUE, prob = assumption.weights)
    
  })
  
  y.true.samp <- unlist(y.true.samp)
  
  y.true.frame <- matrix(data = y.true.samp, ncol = 2, byrow = TRUE)
  
  
  # First Order Estimate 
  train.p1.hat <- c()
  for(y in 0:(N)){
    prop.tmp <- mean(y.true.frame[,1] == y)
    train.p1.hat <- c(train.p1.hat, prop.tmp)
  }
  
  # Second Order Estimate
  p.hat.two.obs <- rep(NA, nrow(y.true.frame))
  two.obs.mapping$two.obs.cat <- 1:nrow(two.obs.mapping)
  for(k in 1:nrow(y.true.frame)){
    idx1 <- two.obs.mapping[,1] == y.true.frame[k,1]
    idx2 <- two.obs.mapping[,2] == y.true.frame[k,2]
    idx <- idx1 & idx2
    p.hat.two.obs[k] <- two.obs.mapping$two.obs.cat[idx]
  }
  
  
  train.p2.hat <- c()
  for(s in 1:((N+1)^2)){
    prop.tmp <- mean(p.hat.two.obs == s)
    train.p2.hat <- c(train.p2.hat, prop.tmp)
  }
  
  
  cat(paste0("Dataset Sample: ", i, " of ", n.sample.sims), end = "\r")
  
  res.list <- lapply(1:nrow(hyper.param.idx), function(x){
    
    j = hyper.param.idx[x,1]
    k = hyper.param.idx[x,2]
    h.tmp <- h.set[j]
    ker.tmp <- ker.set[[k]]
    mu.tmp <- 0
    A.matrix1 <- A.matrix.set.order.1[[j]][[k]]
    A.matrix2 <- A.matrix.set.order.2[[j]][[k]]
    
    
    
    model.estimate1 <- numeric.latent.fit(p.hat = train.p1.hat, 
                                          A.matrix = A.matrix1, 
                                          mu = mu.tmp, 
                                          n.samples = dataset.size,
                                          show.plot = F,
                                          feasibility.test = T)
    
    
    
    
  
    
    # Feasibility may not work for bounded kernels 
    model.estimate2 <-tryCatch({
      numeric.latent.fit(p.hat = train.p2.hat, 
                         A.matrix = A.matrix2, 
                         mu = mu.tmp, 
                         n.samples = dataset.size,
                         show.plot = F,
                         feasibility.test = T)
    },error =function(e){
      placeholder.est()
    }
    )

    
    
    res.tmp <- list("p1" = model.estimate1$p_value, "p2" = model.estimate2$p_value)
    return(res.tmp)
  })
  
  
  
  
  for(q in 1:nrow(hyper.param.idx)){
    j = hyper.param.idx[q,1]
    k = hyper.param.idx[q,2]
    results.array.order.1[i,j,k] <- res.list[[q]]$p1
    results.array.order.2[i,j,k] <- res.list[[q]]$p2
  }
  
}









#saveRDS(results.array.order.1,paste0("data/feasibility_test_order1_sims_200.rds"))
#saveRDS(results.array.order.2,paste0("data/feasibility_test_order2_sims_200.rds"))


results.array.order.1 <- readRDS("data/feasibility_test_order1_sims_200.rds")
results.array.order.2 <- readRDS("data/feasibility_test_order2_sims_200.rds")

n.sample.sims <- nrow(results.array.order.1)
h.set <- c(0.25,0.5,1,2,4,7,10,14,20,25,30)

gauss.frame.1 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))
exp.frame.1 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))
triangle.frame.1 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))
epanechnikov.frame.1 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))

colnames(gauss.frame.1) <- h.set #paste0("h_", h.set)
colnames(exp.frame.1) <- h.set
colnames(triangle.frame.1) <- h.set
colnames(epanechnikov.frame.1) <- h.set


gauss.frame.2 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))
exp.frame.2 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))
triangle.frame.2 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))
epanechnikov.frame.2 <- data.frame(matrix(NA, nrow = n.sample.sims, ncol = length(h.set)))

colnames(gauss.frame.2) <- h.set
colnames(exp.frame.2) <- h.set
colnames(triangle.frame.2) <- h.set
colnames(epanechnikov.frame.2) <- h.set




for(i in 1:n.sample.sims){
  gauss.frame.1[i,] <- results.array.order.1[i,,1]
  exp.frame.1[i,] <- results.array.order.1[i,,2]
  triangle.frame.1[i,] <- results.array.order.1[i,,3]
  epanechnikov.frame.1[i,] <- results.array.order.1[i,,4]
  
  gauss.frame.2[i,] <- results.array.order.2[i,,1]
  exp.frame.2[i,] <- results.array.order.2[i,,2]
  triangle.frame.2[i,] <- results.array.order.2[i,,3]
  epanechnikov.frame.2[i,] <- results.array.order.2[i,,4]
  
}

library(dplyr)
library(tidyr)
gauss.frame.1.long <- gauss.frame.1 %>% gather(bandwidth, p.value)
exp.frame.1.long <- exp.frame.1 %>% gather(bandwidth, p.value)
triangle.frame.1.long <- triangle.frame.1 %>% gather(bandwidth, p.value)
epanechnikov.frame.1.long <- epanechnikov.frame.1 %>% gather(bandwidth, p.value)

gauss.frame.1.long$bandwidth <- factor(gauss.frame.1.long$bandwidth, levels = h.set)
exp.frame.1.long$bandwidth <- factor(exp.frame.1.long$bandwidth, levels = h.set)
triangle.frame.1.long$bandwidth <- factor(triangle.frame.1.long$bandwidth, levels = h.set)
epanechnikov.frame.1.long$bandwidth <- factor(epanechnikov.frame.1.long$bandwidth, levels = h.set)

gauss.frame.2.long <- gauss.frame.2 %>% gather(bandwidth, p.value)
exp.frame.2.long <- exp.frame.2 %>% gather(bandwidth, p.value)
triangle.frame.2.long <- triangle.frame.2 %>% gather(bandwidth, p.value)
epanechnikov.frame.2.long <- epanechnikov.frame.2 %>% gather(bandwidth, p.value)

gauss.frame.2.long$bandwidth <- factor(gauss.frame.2.long$bandwidth, levels = h.set)
exp.frame.2.long$bandwidth <- factor(exp.frame.2.long$bandwidth, levels = h.set)
triangle.frame.2.long$bandwidth <- factor(triangle.frame.2.long$bandwidth, levels = h.set)
epanechnikov.frame.2.long$bandwidth <- factor(epanechnikov.frame.2.long$bandwidth, levels = h.set)


p.gauss.1 <- ggplot(data = gauss.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Gaussian Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

p.exp.1 <- ggplot(data = exp.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Laplace Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

p.triangle.1 <- ggplot(data = triangle.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Triangle Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

p.epanechnikov.1 <- ggplot(data = epanechnikov.frame.1.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Epanechnikov Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  



p.gauss.2 <- ggplot(data = gauss.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Gaussian Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

p.exp.2 <- ggplot(data = exp.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Laplace Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

p.triangle.2 <- ggplot(data = triangle.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Triangle Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  

p.epanechnikov.2 <- ggplot(data = epanechnikov.frame.2.long, aes(x = bandwidth, y = p.value)) + geom_boxplot() + ggtitle("Epanechnikov Kernel Feasibility Tests") + 
  theme(axis.text.x = element_text( size = axis.ticks.size),
        axis.text.y = element_text( size = axis.ticks.size),  
        axis.title.x = element_text( size = axis.size),
        axis.title.y = element_text( size = axis.size),
        title = element_text( size = title.size),
        legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=axis.size), #change legend title font size
        legend.text = element_text(size=axis.ticks.size))  



grid.arrange(arrangeGrob(p.gauss.1, p.exp.1, ncol = 2),                             
             arrangeGrob(p.triangle.1, p.epanechnikov.1, ncol = 2), 
             nrow = 2)   

grid.arrange(arrangeGrob(p.gauss.2, p.exp.2, ncol = 2),                             
             arrangeGrob(p.triangle.2, p.epanechnikov.2, ncol = 2), 
             nrow = 2)   



png(filename = "plots/extrinsic_variability_sim.png",
    width = png.width, height = png.height)



ggarrange(gauss.plot, exp.plot,
          ncol=2, nrow=1, common.legend = TRUE, legend="right") 
# Close the pdf file
dev.off() 



png(filename = "plots/feasibility_sims1_large.png",
    width = png.width, height = png.height)

ggarrange(p.gauss.1, p.exp.1, p.triangle.1, 
          p.epanechnikov.1, nrow = 2, ncol = 2, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 

png(filename = "plots/feasibility_sims1.png",
    width = 2*png.width, height = png.height)

ggarrange(p.gauss.1, p.exp.1, p.epanechnikov.1, nrow = 1, ncol = 3, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 

png(filename = "plots/feasibility_sims2_large.png",
    width = png.width, height = png.height)

ggarrange(p.gauss.2, p.exp.2, p.triangle.2, 
          p.epanechnikov.2, nrow = 2, ncol = 2, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 

png(filename = "plots/feasibility_sims2.png",
    width = 2*png.width, height = png.height)

ggarrange(p.gauss.2, p.exp.2, p.epanechnikov.2, nrow = 1, ncol = 3, 
          common.legend = TRUE, legend="right")   

# Close the pdf file
dev.off() 



