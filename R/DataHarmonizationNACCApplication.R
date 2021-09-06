
# Functions for Data Harmonization
setwd("/Users/Owner/Documents/PhD/Latent Variable Modelling with Application to Data Harmonization")
source("DataHarmonizationFunctions.R")
library(dplyr)
library(gridExtra)

categorize <- function(data){
  tmp <- data %>% mutate(group = ifelse(((sex == 1) & (educ <= 16)), 1, 
                                        ifelse(((sex == 2) & (educ <= 16)), 2, 
                                               ifelse(((sex == 1) & (educ > 16)), 3, 4))))
  tmp <- tmp[, !names(tmp) %in% c('sex','educ')]  
  tmp$group <- factor(tmp$group)
  return(tmp)
}


##### Figure information 
png.width = 1200
png.height = 1000

title.size = 40 #22
axis.ticks.size = 24 #14
axis.size = 30 #18





full.dat <- read.csv("data/investigator_nacc47.csv")
all_tests <- full.dat
w_normal = all_tests$CDRGLOB==0
w_age = (all_tests$NACCAGE>59) & (all_tests$NACCAGE <= 85)

all_tests = all_tests[w_normal&w_age,]
n = nrow(all_tests)

score1 <- "NACCMMSE"
score2 <- "MOCATOTS"

first_score1 <- paste0("first_",score1)
first_score2 <- paste0("first_",score2)

Ny <- 30
Nz <- 30
### select variables (demographic)
all_tests = all_tests %>% mutate(NACCMMSE = ifelse(((NACCMMSE >=0)&(NACCMMSE <= Ny)),NACCMMSE, NA))
all_tests = all_tests %>% mutate(MOCATOTS = ifelse(((MOCATOTS >=0)&(MOCATOTS <= Nz)),MOCATOTS, NA))


all_tests <- all_tests %>% group_by(NACCID) %>% 
  arrange(NACCID,VISITYR,VISITMO,VISITDAY) %>% 
  mutate(row.idx = row_number())
  

col.sub.y <- c("NACCMMSE", "NACCID", "VISITYR","VISITMO","VISITDAY", "row.idx", "NACCAGE", "SEX", "EDUC")
tests.y <- all_tests[,col.sub.y]

col.newnames.y <- c("y", "id", "year","mo","day", "row.idx","age", "sex", "educ")
colnames(tests.y) <- col.newnames.y

tests.y <- tests.y[!is.na(tests.y$y),]

tests.y2 <- tests.y %>% left_join(tests.y[,c("y", "id", "year","mo","day", "row.idx")], by = 'id')
colnames(tests.y2) <- c("y1", "id", "year1","mo1","day1", "row.idx1","age", "sex", "educ", 
                        "y2", "year2","mo2","day2", "row.idx2")

# filtering subsequent visits 
tests.y2 <- tests.y2 %>% filter(row.idx1 + 1 == row.idx2)
tests.y2 <- tests.y2 %>% mutate(date.diff = as.Date(paste0(year2, '/', mo2, '/',day2)) - 
                     as.Date(paste0(year1, '/', mo1, '/',day1)))

# including only differences of at most 500 days
tests.y2 <- tests.y2 %>% filter(date.diff <= 500)
tests.y2 <- tests.y2 %>% group_by(id)  %>%  mutate(row.pair.idx = row_number())
tests.y2 <- tests.y2[tests.y2$row.pair.idx == 1, ]



col.sub.z <- c("MOCATOTS", "NACCID", "VISITYR","VISITMO","VISITDAY", "row.idx", "NACCAGE", "SEX", "EDUC")
tests.z <- all_tests[,col.sub.z]
col.newnames.z <- c("z", "id", "year","mo","day", "row.idx","age", "sex", "educ")
colnames(tests.z) <- col.newnames.z

tests.z <- tests.z[!is.na(tests.z$z),]

tests.z2 <- tests.z %>% left_join(tests.z[,c("z", "id", "year","mo","day", "row.idx")], by = 'id')
colnames(tests.z2) <- c("z1", "id", "year1","mo1","day1", "row.idx1","age", "sex", "educ", 
                        "z2", "year2","mo2","day2", "row.idx2")

# filtering subsequent visits 
tests.z2 <- tests.z2 %>% filter(row.idx1 + 1 == row.idx2)
tests.z2 <- tests.z2 %>% mutate(date.diff = as.Date(paste0(year2, '/', mo2, '/',day2)) - 
                                  as.Date(paste0(year1, '/', mo1, '/',day1)))

# including only differences of at most 500 days
tests.z2 <- tests.z2 %>% filter(date.diff <= 500)
tests.z2 <- tests.z2 %>% group_by(id)  %>%  mutate(row.pair.idx = row_number())
tests.z2 <- tests.z2[tests.z2$row.pair.idx == 1, ]

tests.y2 <- categorize(tests.y2)
tests.z2 <- categorize(tests.z2)

tests.y2.val <- tests.y2[,c('y1','age','y2','group')]
tests.z2.val <- tests.z2[,c('z1','age','z2','group')]

write.csv(tests.y2.val,"data/NACCMMSE_validation.csv", row.names = F)
write.csv(tests.z2.val,"data/MOCATOTS_validation.csv", row.names = F)


# tests.y2 <- as.matrix(as.data.frame(tests.y2[,c('y1','age','y2','group')]))
# tests.z2 <- as.matrix(as.data.frame(tests.z2[,c('z1','age','z2','group')]))
# 
# 
# ny2 <- nrow(tests.y2)
# nz2 <- nrow(tests.z2)
# 
# 
# two.obs.mapping.y <- expand.grid(0:Ny,0:Ny)
# 
# y.pair.obs <- rep(NA, nrow(tests.y2[,c('y1','y2')]))
# two.obs.mapping.y$two.obs.cat <- 1:nrow(two.obs.mapping.y)
# for(k in 1:nrow(tests.y2)){
#   idx1 <- (two.obs.mapping.y[,1] == tests.y2[k,'y1'])
#   idx2 <- (two.obs.mapping.y[,2] == tests.y2[k,'y2'])
#   idx <- idx1 & idx2
#   y.pair.obs[k] <- two.obs.mapping.y$two.obs.cat[idx]
# }
# 
# 
# p.hat.y2 <- c()
# for(i in 1:((Ny+1)^2)){
#   prop.tmp <- mean(y.pair.obs == i)
#   p.hat.y2 <- c(p.hat.y2, prop.tmp)
# }
# 
# 
# 
# 
# two.obs.mapping.z <- expand.grid(0:Nz,0:Nz)
# 
# z.pair.obs <- rep(NA, nrow(tests.z2[,c('z1','z2')]))
# two.obs.mapping.z$two.obs.cat <- 1:nrow(two.obs.mapping.z)
# for(k in 1:nrow(tests.z2)){
#   idx1 <- (two.obs.mapping.z[,1] == tests.z2[k,'z1'])
#   idx2 <- (two.obs.mapping.z[,2] == tests.z2[k,'z2'])
#   idx <- idx1 & idx2
#   z.pair.obs[k] <- two.obs.mapping.z$two.obs.cat[idx]
# }
# 
# 
# p.hat.z2 <- c()
# for(i in 1:((Nz+1)^2)){
#   prop.tmp <- mean(z.pair.obs == i)
#   p.hat.z2 <- c(p.hat.z2, prop.tmp)
# }





y.dat <- read.csv("data/NACCMMSE_train.csv")
z.dat <- read.csv("data/MOCATOTS_train.csv")

y.val <- read.csv("data/NACCMMSE_validation.csv")
z.val <- read.csv("data/MOCATOTS_validation.csv")


colnames(y.dat)[1] <- "y"
colnames(z.dat)[1] <- "z"
Ny <- 30
Nz <- 30


# 
# #### Feasibility test for MMSE 
# 
# obs.set.y <- 0:Ny
# y.summary.dat <- y.dat %>%
#   group_by(y) %>%
#   summarise(Count = n())   
# 
# y.summary.dat <- y.summary.dat[y.summary.dat$y %in% obs.set.y,]
# ny <- sum(y.summary.dat$Count)
# 
# missing.obs <- obs.set.y[!obs.set.y %in%  y.summary.dat$y]
# y.summary.dat <- rbind(y.summary.dat, data.frame(y = missing.obs, Count = rep(0,length(missing.obs))))
# y.summary.dat <- y.summary.dat %>% arrange(y)
# 
# p.hat.y <- y.summary.dat$Count/sum(y.summary.dat$Count)
# 
# pop.y <- numeric.latent.fit(p.hat.y, A1.bin, mu = 0, 
#                             n.samples = ny, show.plot = T, 
#                             feasibility.test = T)
# 
# pop.y$p_value
# 
# 
# 
# 
# #### Feasibility test for MOCA 
# 
# obs.set.z <- 0:Nz
# z.summary.dat <- z.dat %>%
#   group_by(z) %>%
#   summarise(Count = n())   
# 
# z.summary.dat <- z.summary.dat[z.summary.dat$z %in% obs.set.z,]
# nz <- sum(z.summary.dat$Count)
# 
# missing.obs <- obs.set.z[!obs.set.z %in%  z.summary.dat$z]
# z.summary.dat <- rbind(z.summary.dat, data.frame(z = missing.obs, Count = rep(0,length(missing.obs))))
# z.summary.dat <- z.summary.dat %>% arrange(z)
# 
# p.hat.z <- z.summary.dat$Count/sum(z.summary.dat$Count)
# 
# pop.z <- numeric.latent.fit(p.hat.z, A1.bin, mu = 0, 
#                             n.samples = nz, show.plot = T, 
#                             feasibility.test = T)
# 
# pop.z$p_value
# 
# 
# 



## 2nd order population feasibility 
# 
# A2.bin <- A.matrix.binomial.two.obs(R_bins = 1000, N = 30, two.obs.mapping = two.obs.mapping.y, numeric.points = 400)
# 
# 
# pop.y2 <- numeric.latent.fit(p.hat.y2, A2.bin, mu = 0, 
#                             n.samples = ny2, show.plot = T, 
#                             feasibility.test = T)
# 
# pop.y2$p_value
# 
# 
# 
# pop.z2 <- numeric.latent.fit(p.hat.z2, A2.bin, mu = 0, 
#                              n.samples = nz2, show.plot = T, 
#                              feasibility.test = T)
# 
# 
# pop.z2$p_value



age.range <- 60:85
group.range <- 1:4
x.design <- data.frame(age = rep(age.range, times = length(group.range)), 
                       group = rep(group.range, each = length(age.range)))

p.hat.y.set <- cond.dist.est(x.design = x.design, train.data = y.dat, outcome = "y", N = 30, age.window = 3)


wide.p.hat <- cbind(x.design,p.hat.y.set)
colnames(wide.p.hat) <- c("age", "group", 0:Ny)

wide.p.hat$model.idx = 1:nrow(wide.p.hat)

#### matching rows of the validation set to the observed data conditional estimate 

y.val <- y.val %>% left_join(wide.p.hat[,c('age', 'group','model.idx')], by = c('age', 'group')) 



######## Intrinsic Variability Matching ######### 

h.set<- exp(seq(log(0.8),log(15), length.out = 15))
ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 2)))



results.array <- array(NA, dim = c(length(h.set), length(ker.set), length(mu.set)))
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))


# remove unregularized compact support pairs (becomes uncomputable)
idx1 <- hyper.param.idx[,1] <= 1.0  
idx2 <- hyper.param.idx[,2] %in% c(3,4) 

# Also remove compact kernels with h < 1

idx3 <- hyper.param.idx[,1] < 1 


hyper.param.idx <- hyper.param.idx[!(idx1 & idx2 ), ]
# Pre computing A matrices
A.matrix.set <- list()
for(j in 1:length(h.set)){
  h.tmp <- h.set[j]
  A.matrix.set.tmp <- list()
  for(k in 1:length(ker.set)){
    
    ker.tmp <- ker.set[[k]]
    
    A.mat <- A.matrix.compute(R_bins = 1000, N = Ny, ker = ker.tmp, h = h.tmp, numeric.points = 400)
    A.matrix.set.tmp[[k]] <- A.mat 
    
  }
  A.matrix.set[[j]] <- A.matrix.set.tmp
}


##########################################
hist.cols <- as.character(0:Ny)




train.p.hat <- c()
for(y in 0:(Ny)){
  prop.tmp <- mean(y.val[,1] == y)
  train.p.hat <- c(train.p.hat, prop.tmp)
}


res.list <- lapply(1:nrow(hyper.param.idx), function(x){
  j = hyper.param.idx[x,1]
  k = hyper.param.idx[x,2]
  l = hyper.param.idx[x,3]
  h.tmp <- h.set[j]
  ker.tmp <- ker.set[[k]]
  mu.tmp <- mu.set[l]
  A.matrix <- A.matrix.set[[j]][[k]]
  

  model.estimate <- numeric.latent.fit(p.hat = train.p.hat, 
                                       A.matrix = A.matrix, 
                                       mu = mu.tmp, show.plot = F)
  
  
  latent.mix.list <- list()
  model.observed.list <- list()
  
  
  for (m in 1:nrow(y.val)) {
    latent.mix.list[[m]] <- model.estimate$latent
    model.observed.list[[m]] <- model.estimate$observed
  }
  

  
  res.tmp <- intrinsic.variability(y.true.frame = y.val[,c(1,3)], latent.mix.list = latent.mix.list, 
                                   model.observed.list = model.observed.list, n.samp = 5, N = Ny,
                                   ker = ker.tmp, h = h.tmp, parallel = F, show.plot = F)
  
  return(res.tmp)
})





res.set <- unlist(res.list)
for(q in 1:nrow(hyper.param.idx)){
  j = hyper.param.idx[q,1]
  k = hyper.param.idx[q,2]
  l = hyper.param.idx[q,3]
  results.array[j,k,l] <- res.set[q]
}



#saveRDS(results.array,"data/intrinsic_variability_MMSE.rds")
results.array <- readRDS("data/intrinsic_variability_MMSE.rds")
############# Binomial Model 
R_bins = 1000
A.binom <- A.matrix.binomial(R_bins, N = Ny, numeric.points = 400)


latent.mix.list <- list()
model.observed.list <- list()
mu.tmp <- 0.001

model.estimate <- numeric.latent.fit(p.hat = train.p.hat, 
                                     A.matrix = A.binom, 
                                     mu = mu.tmp, show.plot = F)


latent.mix.list <- list()
model.observed.list <- list()


for (m in 1:nrow(y.val)) {
  latent.mix.list[[m]] <- model.estimate$latent
  model.observed.list[[m]] <- model.estimate$observed
}





binom.tot.var <- intrinsic.variability.binom(y.true.frame = y.val[,c(1,3)], latent.mix.list = latent.mix.list, 
                                             model.observed.list = model.observed.list, n.samp = 20, N = Ny,
                                             parallel = F, show.plot = T)




##########################################







mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(mu.set))

gauss.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,1,]))
exp.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,2,]))
triangle.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,3,]))
epanechnikov.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,4,]))


# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Gaussian Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                        axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                        axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                        axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                        title = element_text( size = title.size),
                                                                                                                                                                                                                        legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                        legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                        legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                        legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                        legend.text = element_text(size=axis.ticks.size))
exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Laplace Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                   axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                   axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                   axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                   title = element_text( size = title.size),
                                                                                                                                                                                                                   legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                   legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                   legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                   legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                   legend.text = element_text(size=axis.ticks.size))
triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Triangle Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                               axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                               axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                               axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                               title = element_text( size = title.size),
                                                                                                                                                                                                                               legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                               legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                               legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                               legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                               legend.text = element_text(size=axis.ticks.size))
epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Epanechnikov Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                           axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                           axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                           axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                           title = element_text( size = title.size),
                                                                                                                                                                                                                                           legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                           legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                           legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                           legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                           legend.text = element_text(size=axis.ticks.size))


gauss.plot
exp.plot
triangle.plot
epanechnikov.plot

grid.arrange(arrangeGrob(gauss.plot, exp.plot, ncol = 2),                             
             arrangeGrob(triangle.plot, epanechnikov.plot, ncol = 2), 
             nrow = 2)   


png(filename = "plots/Intrinsic_Variability_MMSE.png",
    width = png.width, height = png.height)

ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 


min(exp.frame$total.variation)
binom.tot.var

exp.frame$h[exp.frame$total.variation == min(exp.frame$total.variation)]

#### MMSE chooses the binomial as the best model 
results.array.y <- readRDS("data/intrinsic_variability_MMSE.rds")










#### Mu Selection 





############# Laplace Model 
mu.set <- c(0,exp(seq(log(0.0001),log(0.1), length.out = 20)))
results.array.mu.selection <- rep(0, length(mu.set))

R_bins = 1000


A.model.y <- A.matrix.compute(R_bins, N = Ny, 
                              ker = ker.y, h = hy,
                              numeric.points = 400)

A.two.sample.model.y <- A.two.obs.tensor.compute(R_bins = 1000, 
                                                 N = Ny, ker = ker.y, h = hy,
                                                 numeric.points = 400)


latent.mix.list.short <- list()
model.observed.list.short <- list()

latent.mix.list <- list()
model.observed.list <- list()
mu.tmp <- 0.001





res.list <- lapply(1:length(mu.set), function(x){

  mu.tmp <- mu.set[x]
  
  for (m in 1:nrow(wide.p.hat)) {
    
    train.p.hat.x <- as.numeric(wide.p.hat[m,hist.cols])
    model.estimate <- numeric.latent.fit(p.hat = train.p.hat.x, 
                                         A.matrix = A.model.y, 
                                         mu = mu.tmp, show.plot = F)
    
    latent.mix.list.short[[m]] <- model.estimate$latent
    model.observed.list.short[[m]] <- model.estimate$observed
    cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat)), end="\r")
  }
  cat(end="\n")
  for(m in 1:nrow(y.val)){
    k = y.val$model.idx[m]
    latent.mix.list[[m]] <-latent.mix.list.short[[k]]
    model.observed.list[[m]] <- model.observed.list.short[[k]]
    
  }
  
 
  
  
  log.like <- two.samp.log.likelihood(y.val[,c(1,3)], 
                                      A.two.sample.model.y, 
                                      latent.mix.list)
  
  
  res.tmp <- -log.like
  return(res.tmp)
})

res.set <- unlist(res.list)

for(q in 1:length(mu.set)){
  results.array.mu.selection[q] <- res.set[q]
}


#saveRDS(results.array.mu.selection,paste0("data/mu_selection_MMSE.rds"))
results.array.mu.selection <- readRDS("data/mu_selection_MMSE.rds")

mu.frame <- data.frame(log.mu = log(mu.set + 0.001), neg.likelihood = results.array.mu.selection)
mu.plot <- ggplot(data = mu.frame, aes(x = log.mu, y = neg.likelihood)) + geom_line() + ylab("Negative Two Sample Log-Likelihood")+ ggtitle("Regularization With Selected Binomial Model") + xlab("log(\u00b5 + 0.001)") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                  axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                  axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                  axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                  title = element_text( size = title.size),
                                                                                                                                                                                                                                  legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                  legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                  legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                  legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                  legend.text = element_text(size=axis.ticks.size))

png(filename = "plots/mu_selection_MMSE.png",
    width = png.width, height = png.height)

mu.plot

# Close the pdf file
dev.off() 


print(paste0("Optimal Mu Value: ", mu.set[which.min(results.array.mu.selection)]))

opt.mu.y <- mu.set[which.min(results.array.mu.selection)]

#### Repeat for Z






age.range <- 60:85
group.range <- 1:4
x.design <- data.frame(age = rep(age.range, times = length(group.range)), 
                       group = rep(group.range, each = length(age.range)))

p.hat.z.set <- cond.dist.est(x.design = x.design, train.data = z.dat, outcome = "z", N = Nz, age.window = 3)


wide.p.hat <- cbind(x.design,p.hat.z.set)
colnames(wide.p.hat) <- c("age", "group", 0:Nz)

wide.p.hat$model.idx = 1:nrow(wide.p.hat)

#### matching rows of the validation set to the observed data conditional estimate 

z.val <- z.val %>% left_join(wide.p.hat[,c('age', 'group','model.idx')], by = c('age', 'group')) 




######## Intrinsic Variability Matching ######### 

h.set<- exp(seq(log(0.8),log(15), length.out = 15))
ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.001),log(0.1), length.out = 2)))



results.array <- array(NA, dim = c(length(h.set), length(ker.set), length(mu.set)))
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))


# remove unregularized compact support pairs (becomes uncomputable)
idx1 <- hyper.param.idx[,1] <= 1.0  
idx2 <- hyper.param.idx[,2] %in% c(3,4) 

# Also remove compact kernels with h < 1

idx3 <- hyper.param.idx[,1] < 1 


hyper.param.idx <- hyper.param.idx[!(idx1 & idx2 ), ]
# Pre computing A matrices
A.matrix.set <- list()
for(j in 1:length(h.set)){
  h.tmp <- h.set[j]
  A.matrix.set.tmp <- list()
  for(k in 1:length(ker.set)){
    
    ker.tmp <- ker.set[[k]]
    
    A.mat <- A.matrix.compute(R_bins = 1000, N = Ny, ker = ker.tmp, h = h.tmp, numeric.points = 400)
    A.matrix.set.tmp[[k]] <- A.mat 
    
  }
  A.matrix.set[[j]] <- A.matrix.set.tmp
}


##########################################
hist.cols <- as.character(0:Nz)




train.p.hat <- c()
for(z in 0:(Nz)){
  prop.tmp <- mean(z.val[,1] == z)
  train.p.hat <- c(train.p.hat, prop.tmp)
}


res.list <- lapply(1:nrow(hyper.param.idx), function(x){
  j = hyper.param.idx[x,1]
  k = hyper.param.idx[x,2]
  l = hyper.param.idx[x,3]
  h.tmp <- h.set[j]
  ker.tmp <- ker.set[[k]]
  mu.tmp <- mu.set[l]
  A.matrix <- A.matrix.set[[j]][[k]]
  
  
  model.estimate <- numeric.latent.fit(p.hat = train.p.hat, 
                                       A.matrix = A.matrix, 
                                       mu = mu.tmp, show.plot = F)
  
  
  latent.mix.list <- list()
  model.observed.list <- list()
  
  
  for (m in 1:nrow(z.val)) {
    latent.mix.list[[m]] <- model.estimate$latent
    model.observed.list[[m]] <- model.estimate$observed
  }
  
  
  
  res.tmp <- intrinsic.variability(y.true.frame = z.val[,c(1,3)], latent.mix.list = latent.mix.list, 
                                   model.observed.list = model.observed.list, n.samp = 5, N = Nz,
                                   ker = ker.tmp, h = h.tmp, parallel = F, show.plot = F)
  
  return(res.tmp)
})





res.set <- unlist(res.list)
for(q in 1:nrow(hyper.param.idx)){
  j = hyper.param.idx[q,1]
  k = hyper.param.idx[q,2]
  l = hyper.param.idx[q,3]
  results.array[j,k,l] <- res.set[q]
}


#saveRDS(results.array,"data/intrinsic_variability_MOCA.rds")
results.array <- readRDS("data/intrinsic_variability_MOCA.rds")


############# Binomial Model 
R_bins = 1000
A.binom <- A.matrix.binomial(R_bins, N = Nz, numeric.points = 400)


latent.mix.list <- list()
model.observed.list <- list()
mu.tmp <- 0.001

model.estimate <- numeric.latent.fit(p.hat = train.p.hat, 
                                     A.matrix = A.binom, 
                                     mu = mu.tmp, show.plot = F)


latent.mix.list <- list()
model.observed.list <- list()


for (m in 1:nrow(z.val)) {
  latent.mix.list[[m]] <- model.estimate$latent
  model.observed.list[[m]] <- model.estimate$observed
}





binom.tot.var <- intrinsic.variability.binom(y.true.frame = z.val[,c(1,3)], latent.mix.list = latent.mix.list, 
                                             model.observed.list = model.observed.list, n.samp = 2, N = Nz,
                                             parallel = F, show.plot = T)




##########################################







mus <- round(rep(mu.set, each = length(h.set)), 3)
hs <- rep(h.set, times = length(mu.set))

gauss.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,1,]))
exp.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,2,]))
triangle.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,3,]))
epanechnikov.frame <- data.frame(h = hs, mu = mus, total.variation = as.numeric(results.array[,4,]))


# argmin is the gaussian kernel with mu = 0.01

gauss.frame$mu <- factor(gauss.frame$mu)
exp.frame$mu <- factor(exp.frame$mu)
triangle.frame$mu <- factor(triangle.frame$mu)
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu)







gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Gaussian Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                        axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                        axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                        axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                        title = element_text( size = title.size),
                                                                                                                                                                                                                        legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                        legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                        legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                        legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                        legend.text = element_text(size=axis.ticks.size))
exp.plot <- ggplot(data = exp.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line() + ggtitle("Laplace Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                   axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                   axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                   axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                   title = element_text( size = title.size),
                                                                                                                                                                                                                   legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                   legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                   legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                   legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                   legend.text = element_text(size=axis.ticks.size))
triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Triangle Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                               axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                               axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                               axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                               title = element_text( size = title.size),
                                                                                                                                                                                                                               legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                               legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                               legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                               legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                               legend.text = element_text(size=axis.ticks.size))
epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = total.variation, color = mu, group = mu)) + geom_line()  + ggtitle("Epanechnikov Kernel") + geom_hline(yintercept = binom.tot.var) + ylab("Total Variation") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                           axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                           axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                           axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                           title = element_text( size = title.size),
                                                                                                                                                                                                                                           legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                           legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                           legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                           legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                           legend.text = element_text(size=axis.ticks.size))


gauss.plot
exp.plot
triangle.plot
epanechnikov.plot

grid.arrange(arrangeGrob(gauss.plot, exp.plot, ncol = 2),                             
             arrangeGrob(triangle.plot, epanechnikov.plot, ncol = 2), 
             nrow = 2)   



png(filename = "plots/Intrinsic_Variability_MOCA.png",
    width = png.width, height = png.height)

ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 

min(exp.frame$total.variation)
min(gauss.frame$total.variation)
binom.tot.var

exp.frame$h[exp.frame$total.variation == min(exp.frame$total.variation)]











#### Mu Selection 





############# Laplace Model 
mu.set <- c(0,exp(seq(log(0.0001),log(0.1), length.out = 20)))
results.array.mu.selection <- rep(0, length(mu.set))

R_bins = 1000
A.laplace <- A.matrix.compute(R_bins, N = Nz, 
                              ker = exponential_kernel, h = 1.5,
                              numeric.points = 400)


A.two.sample.laplace <- A.two.obs.tensor.compute(R_bins = 1000, 
                                                N = Nz, ker = exponential_kernel, h = 1.5,
                                                numeric.points = 400)


latent.mix.list.short <- list()
model.observed.list.short <- list()

latent.mix.list <- list()
model.observed.list <- list()






res.list <- lapply(1:length(mu.set), function(x){
  
  mu.tmp <- mu.set[x]
  
  for(m in 1:nrow(wide.p.hat)) {
    
    train.p.hat.x <- as.numeric(wide.p.hat[m,hist.cols])
    model.estimate <- numeric.latent.fit(p.hat = train.p.hat.x, 
                                         A.matrix = A.laplace, 
                                         mu = mu.tmp, show.plot = F)
    
    latent.mix.list.short[[m]] <- model.estimate$latent
    model.observed.list.short[[m]] <- model.estimate$observed
    cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat)), end="\r")
  }
  cat(end="\n")
  for(m in 1:nrow(z.val)){
    k = z.val$model.idx[m]
    latent.mix.list[[m]] <-latent.mix.list.short[[k]]
    model.observed.list[[m]] <- model.observed.list.short[[k]]
    
  }
  
  
  log.like <- two.samp.log.likelihood(z.val[,c(1,3)], 
                                      A.two.sample.laplace, 
                                      latent.mix.list)
  
  
  res.tmp <- -log.like
  return(res.tmp)
})

res.set <- unlist(res.list)

for(q in 1:length(mu.set)){
  results.array.mu.selection[q] <- res.set[q]
}


# saveRDS(results.array.mu.selection,paste0("data/mu_selection_MOCA.rds"))
results.array.mu.selection <- readRDS("data/mu_selection_MOCA.rds")

mu.frame <- data.frame(log.mu = log(mu.set + 0.001), neg.likelihood = results.array.mu.selection)
mu.plot <- ggplot(data = mu.frame, aes(x = log.mu, y = neg.likelihood)) + geom_line() + ylab("Negative Two Sample Log-Likelihood")+ ggtitle("Regularization With Selected Laplace Model") + xlab("log(\u00b5 + 0.001)") + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                 axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                 axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                 axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                 title = element_text( size = title.size),
                                                                                                                                                                                                                                 legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                 legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                 legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                 legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                 legend.text = element_text(size=axis.ticks.size))
mu.plot


png(filename = "plots/mu_selection_MOCA.png",
    width = png.width, height = png.height)

mu.plot
# Close the pdf file
dev.off() 

print(paste0("Optimal Mu Value: ", mu.set[which.min(results.array.mu.selection)]))

opt.mu.z <- mu.set[which.min(results.array.mu.selection)]









##### Conversion Cross Entropy: 

# Once again we fix the learned MMSE model and compute the test cross entropy 

cw.test <- read.csv(file = "data/NACCMMSE_to_MOCATOTS_test.csv")
nrow(cw.test)







############################################################
############################################################

# We now test the conversion error 
opt.mu.y <- 0.002
opt.mu.z <- 0.0008

mu.y <- opt.mu.y # 0.002
mu.z <- opt.mu.z # 0.0008


# Assume we have already computed an optimal h ker and mu from the first model 
Ny = 30
Nz = 30

num.bins <- 300


A.model.y <- A.matrix.binomial(R_bins = num.bins, N = Ny, 
                               numeric.points = 400)

h.set<- c(0.35,1.0,1.5,2.0,2.5,4.0,6.0,8.0,12.0) 

ker.set <- list(gaussian_kernel, exponential_kernel, triangle_kernel, epanechnikov_kernel)
mu.set <- c(0,exp(seq(log(0.0001),log(mu.z), length.out = 3)), exp(seq(log(mu.z),log(.3), length.out = 5)))
mu.set <- unique(mu.set)

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




results.array <- array(NA, dim = c(length(h.set), length(ker.set), length(mu.set)))
results.array.ml <- array(NA, dim = c(length(h.set), length(ker.set)))
results.array.para <- results.array.ml
hyper.param.idx <- expand.grid(1:length(h.set),1:length(ker.set), 1:length(mu.set))
hyper.param.idx.ml <- expand.grid(1:length(h.set),1:length(ker.set))


# compact support kernels must have h >= 1

idx1 <- hyper.param.idx[,1]  %in% which(h.set <= 1)
idx2 <- hyper.param.idx[,2]  %in% c(3,4)

hyper.param.idx <- hyper.param.idx[!(idx1 & idx2),]


idx1 <- hyper.param.idx.ml[,1]  %in% which(h.set <= 1)
idx2 <- hyper.param.idx.ml[,2]  %in% c(3,4)

hyper.param.idx.ml <- hyper.param.idx.ml[!(idx1 & idx2),]


### population test set 






age.range <- 60:85
group.range <- 1:4
x.design <- data.frame(age = rep(age.range, times = length(group.range)), 
                       group = rep(group.range, each = length(age.range)))

p.hat.y.set <- cond.dist.est(x.design = x.design, train.data = y.dat, outcome = "y", N = Ny, age.window = 3)
p.hat.z.set <- cond.dist.est(x.design = x.design, train.data = z.dat, outcome = "z", N = Nz, age.window = 3)


wide.p.hat.y <- cbind(x.design,p.hat.y.set)
colnames(wide.p.hat.y) <- c("age", "group", 0:Ny)

wide.p.hat.z <- cbind(x.design,p.hat.z.set)
colnames(wide.p.hat.z) <- c("age", "group", 0:Nz)



wide.p.hat.y$model.idx = 1:nrow(wide.p.hat.y)
wide.p.hat.z$model.idx = 1:nrow(wide.p.hat.z)

#### matching rows of the validation set to the observed data conditional estimate 

cw.test <- cw.test %>% left_join(wide.p.hat.y[,c('age', 'group','model.idx')], by = c('age', 'group')) 


hist.cols.y <- as.character(0:Ny)
hist.cols.z <- as.character(0:Nz)

latent.mix.list.short.y <- list()
model.observed.list.short.y <- list()
latent.mix.list.short.y.ml <- list()
model.observed.list.short.y.ml <- list()


for(m in 1:nrow(wide.p.hat.y)){
  train.p.hat.y <- as.numeric(wide.p.hat.y[m,hist.cols.y])
  model.estimate.y <- numeric.latent.fit(p.hat = train.p.hat.y, 
                                       A.matrix = A.model.y, 
                                       mu = opt.mu.y, show.plot = F)
  
  latent.mix.list.short.y[[m]] <- model.estimate.y$latent
  model.observed.list.short.y[[m]] <- model.estimate.y$observed
  
  model.estimate.y.ml <- numeric.latent.fit(p.hat = train.p.hat.y, 
                                            A.matrix = A.model.y, 
                                            mu = 0, show.plot = F)
  
  latent.mix.list.short.y.ml[[m]] <- model.estimate.y.ml$latent
  model.observed.list.short.y.ml[[m]] <- model.estimate.y.ml$observed
  
  
  cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat.y)), end="\r")
}

cat(end="\n")

latent.mix.list.y <- list()
model.observed.list.y <- list()

latent.mix.list.y.ml <- list()
model.observed.list.y.ml <- list()

for(m in 1:nrow(cw.test)){
  k = cw.test$model.idx[m]
  latent.mix.list.y[[m]] <-latent.mix.list.short.y[[k]]
  model.observed.list.y[[m]] <- model.observed.list.short.y[[k]]
  
  latent.mix.list.y.ml[[m]] <- latent.mix.list.short.y.ml[[k]]
  model.observed.list.y.ml[[m]] <- model.observed.list.short.y.ml[[k]]
}



### Fixed Model Estimate of y 




res <- sapply(1:nrow(hyper.param.idx), function(x){
  j = hyper.param.idx[x,1]
  k = hyper.param.idx[x,2]
  l = hyper.param.idx[x,3]
  h.tmp <- h.set[j]
  ker.tmp <- ker.set[[k]]
  mu.tmp <- mu.set[l]
  A.matrix <- A.matrix.set[[j]][[k]]

  
  ## fitting Latent Models 
  
  latent.mix.list.short.z <- list()
  model.observed.list.short.z <- list()
  
  
  for(m in 1:nrow(wide.p.hat.z)){
    train.p.hat.z <- as.numeric(wide.p.hat.z[m,hist.cols.z])
    model.estimate.z <- numeric.latent.fit(p.hat = train.p.hat.z, 
                                           A.matrix = A.matrix, 
                                           mu = mu.tmp, show.plot = F)
    
    latent.mix.list.short.z[[m]] <- model.estimate.z$latent
    model.observed.list.short.z[[m]] <- model.estimate.z$observed
    
    
    cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat.z)), end="\r")
  }
  
  cat(end="\n")
  
  latent.mix.list.z <- list()
  model.observed.list.z <- list()
  
  
  for(m in 1:nrow(cw.test)){
    k = cw.test$model.idx[m]
    latent.mix.list.z[[m]] <-latent.mix.list.short.z[[k]]
    model.observed.list.z[[m]] <- model.observed.list.short.z[[k]]
    
  }
  
  
  
  res.tmp <- score.conversion(test.pairs = cw.test[,c(1,2)], 
                              latent.mix.list.y = latent.mix.list.y,
                              model.observed.list.y = model.observed.list.y,
                              latent.mix.list.z = latent.mix.list.z, 
                              model.observed.list.z = model.observed.list.z, 
                              Ny = Ny, 
                              Nz = Nz, ker.z = ker.tmp, hz = h.tmp, 
                              n.samp = 4, method = "CrossEntropy", binomial.y = T)
  
  
  
  
  cat(paste0("Param Set: ", x, "/", nrow(hyper.param.idx)), end = "\r")
  return(res.tmp)
})
cat( end = "\n")
res.set <- res


for(q in 1:nrow(hyper.param.idx)){
  j = hyper.param.idx[q,1]
  k = hyper.param.idx[q,2]
  l = hyper.param.idx[q,3]
  results.array[j,k,l] <- res.set[q]
}


cat(paste0("ML conversion"), end = "\n")

##### ML conversion 

res <- sapply(1:nrow(hyper.param.idx.ml), function(x){
  j = hyper.param.idx.ml[x,1]
  k = hyper.param.idx.ml[x,2]
  h.tmp <- h.set[j]
  ker.tmp <- ker.set[[k]]
  
  A.matrix <- A.matrix.set[[j]][[k]]
  
  
  latent.mix.list.short.z.ml <- list()
  model.observed.list.short.z.ml <- list()
  
  
  for(m in 1:nrow(wide.p.hat.z)){
    train.p.hat.z <- as.numeric(wide.p.hat.z[m,hist.cols.z])
    model.estimate.z.ml <- numeric.latent.fit(p.hat = train.p.hat.z, 
                                              A.matrix = A.matrix, 
                                              mu = 0, show.plot = F)
    
    latent.mix.list.short.z.ml[[m]] <- model.estimate.z.ml$latent
    model.observed.list.short.z.ml[[m]] <- model.estimate.z.ml$observed
    
    
    cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat.z)), end="\r")
  }
  
  cat(end="\n")
  
  latent.mix.list.z.ml <- list()
  model.observed.list.z.ml <- list()
  
  
  for(m in 1:nrow(cw.test)){
    k = cw.test$model.idx[m]
    latent.mix.list.z.ml[[m]] <-latent.mix.list.short.z.ml[[k]]
    model.observed.list.z.ml[[m]] <- model.observed.list.short.z.ml[[k]]
    
  }
  
  res.tmp <- score.conversion(test.pairs = cw.test[,c(1,2)], 
                              latent.mix.list.y = latent.mix.list.y.ml,
                              model.observed.list.y = model.observed.list.y.ml,
                              latent.mix.list.z = latent.mix.list.z.ml, 
                              model.observed.list.z = model.observed.list.z.ml, 
                              Ny = Ny, 
                              Nz = Nz, ker.z = ker.tmp, hz = h.tmp, 
                              n.samp = 4, method = "CrossEntropy", binomial.y = T)
  
  
  
  
  
  cat(paste0("Param Set: ", x, "/", nrow(hyper.param.idx.ml)), end = "\r")
  return(res.tmp)
})
cat( end = "\n")
res.set <- res
for(q in 1:nrow(hyper.param.idx.ml)){
  j = hyper.param.idx.ml[q,1]
  k = hyper.param.idx.ml[q,2]
  results.array.ml[j,k] <- res.set[q]
}



cat(paste0("Parametric conversion"), end = "\n")

##### Parametric conversion 
latent.variational.moments.y <- renormalized_r_function(N, is.binomial = T)
params.y <- fit_logitnorm_model(y.dat.wide$y, X = y.dat.wide[,2:6], latent.variational.moments.y)

res <- sapply(1:nrow(hyper.param.idx.ml), function(x){
  j = hyper.param.idx.ml[x,1]
  k = hyper.param.idx.ml[x,2]
  h.tmp <- h.set[j]
  ker.tmp <- ker.set[[k]]
  
  
  latent.variational.moments.z <- renormalized_r_function(N, ker = ker.tmp, h = h.tmp)
  
  params.z <- fit_logitnorm_model(z.dat.wide$z, X = z.dat.wide[,2:6], latent.variational.moments.z)
  
  
  conv.probs <- rep(NA, nrow(cw.test))
  for(m in 1:nrow(cw.test)){
    y = cw.test$NACCMMSE[m]
    z = cw.test$MOCATOTS[m]
    x.row = cw.test.wide[m,c("age", "group.1", "group.2","group.3", "group.4")]
    cond.prob <- conversion.parametric(y = y, z = z, x = x.row, params.y = params.y, params.z = params.z,
                                       Ny = Ny, Nz = Nz, ker.z = ker.tmp, hz = h.tmp, binomial.y = T, R_bins = 1000)
    conv.probs[m] <- cond.prob
    #cat(paste0("Test Sample ", m , "/", nrow(cw.test)), end = "\r")
  }
  #cat( end = "\n")
  
  
  res.tmp <- -sum(log(conv.probs))
  
  
  
  cat(paste0("Param Set: ", x, "/", nrow(hyper.param.idx.ml)), end = "\r")
  #cat( end = "\n")
  return(res.tmp)
})
cat( end = "\n")
res.set <- res
for(q in 1:nrow(hyper.param.idx.ml)){
  j = hyper.param.idx.ml[q,1]
  k = hyper.param.idx.ml[q,2]
  results.array.para[j,k] <- res.set[q]
}



### Binomial Conversion 
A.matrix <- A.matrix.binomial(R_bins = num.bins, 
                              N = Nz,numeric.points = 400)
binom.results <- rep(NA, length(mu.set))
for(j in 1:length(mu.set)){
  mu.tmp <- mu.set[j]
  
  
  latent.mix.list.short.z <- list()
  model.observed.list.short.z <- list()
  
  
  for(m in 1:nrow(wide.p.hat.z)){
    train.p.hat.z <- as.numeric(wide.p.hat.z[m,hist.cols.z])
    model.estimate.z <- numeric.latent.fit(p.hat = train.p.hat.z, 
                                           A.matrix = A.matrix, 
                                           mu = mu.tmp, show.plot = F)
    
    latent.mix.list.short.z[[m]] <- model.estimate.z$latent
    model.observed.list.short.z[[m]] <- model.estimate.z$observed
    
    
    cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat.z)), end="\r")
  }
  
  cat(end="\n")
  
  latent.mix.list.z <- list()
  model.observed.list.z <- list()
  
  
  for(m in 1:nrow(cw.test)){
    k = cw.test$model.idx[m]
    latent.mix.list.z[[m]] <-latent.mix.list.short.z[[k]]
    model.observed.list.z[[m]] <- model.observed.list.short.z[[k]]
    
  }
  
  res.binom <- score.conversion(test.pairs = cw.test[,c(1,2)], 
                                latent.mix.list.y = latent.mix.list.y,
                                model.observed.list.y = model.observed.list.y,
                                latent.mix.list.z = latent.mix.list.z, 
                                model.observed.list.z = model.observed.list.z, 
                                Ny = Ny, 
                                Nz = Nz, 
                                n.samp = 4, method = "CrossEntropy", 
                                binomial.y = T, binomial.z = T)
  
  binom.results[j] <- res.binom
  cat(paste0("Mu: ", j,"/",length(mu.set)), end="\n")
}



# naive method: 

muy <- mean(y.dat$y)
sdy <- sd(y.dat$y)
muz <- mean(z.dat$z)
sdz <- sd(z.dat$z)

naive.conversion.prob <- function(y,muy,sdy,muz,sdz, Nz, n.samp = 100000){
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


ce = 0
for(i in 1:nrow(cw.test)){
  z.prob = naive.conversion.prob(y = cw.test[i,1], muy = muy, sdy = sdy, 
                                 muz = muz, sdz = sdz, Nz = Nz, n.samp = 100000)
  ce = ce -log(z.prob[cw.test[i,2] + 1])
  cat(paste0("Row ", i, "/",nrow(cw.test)), end = "\r")
}








#saveRDS(results.array, "data/conversion_ce_MMSE_MOCA.rds")
#saveRDS(results.array.ml, "data/conversion_ce_ml_MMSE_MOCA.rds")
#saveRDS(results.array.para, "data/conversion_ce_para_MMSE_MOCA.rds")

results.array <- readRDS("data/conversion_ce_MMSE_MOCA.rds")
results.array.ml <- readRDS("data/conversion_ce_ml_MMSE_MOCA.rds")
results.array.para <- readRDS("data/conversion_ce_para_MMSE_MOCA.rds")

mus <- round(rep(mu.set, each = length(h.set)), 5)
mus <- as.character(mus)
mus <- c(mus, rep("Binomial", length(h.set)), rep("Parametric", length(h.set)), rep("Unregularized", length(h.set)), rep("Z-Score Matching", length(h.set)))
hs <- rep(h.set, times = length(mu.set) + 4)

gauss.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(results.array[,1,]), rep(min(binom.results), length(h.set)),as.vector(results.array.para[,1]) ,as.vector(results.array.ml[,1]), rep(ce, length(h.set)) ) )
exp.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(results.array[,2,]), rep(min(binom.results), length(h.set)),as.vector(results.array.para[,2]) ,as.vector(results.array.ml[,2]), rep(ce, length(h.set)) ) )
triangle.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(results.array[,3,]), rep(min(binom.results), length(h.set)),as.vector(results.array.para[,3]) ,as.vector(results.array.ml[,3]), rep(ce, length(h.set)) ) )
epanechnikov.frame <- data.frame(h = hs, mu = mus, cross.entropy = c(as.vector(results.array[,4,]), rep(min(binom.results), length(h.set)),as.vector(results.array.para[,4]) ,as.vector(results.array.ml[,4]), rep(ce, length(h.set)) ) )




gauss.frame$mu <- factor(gauss.frame$mu, levels = as.character(c(round(mu.set, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))
exp.frame$mu <- factor(exp.frame$mu, levels = as.character(c(round(mu.set, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))
triangle.frame$mu <- factor(triangle.frame$mu, levels = as.character(c(round(mu.set, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))
epanechnikov.frame$mu <- factor(epanechnikov.frame$mu, levels = as.character(c(round(mu.set, 5), "Binomial", "Parametric", "Unregularized", "Z-Score Matching")))


mus.small <- round(rep(mu.set, each = length(h.set)), 5)
mus.small <- as.character(mus.small)
mus.small <- c(mus.small, rep("Parametric", length(h.set)), rep("Unregularized", length(h.set)), rep("Z-Score Matching", length(h.set)))
hs.small <- rep(h.set, times = length(mu.set) + 3)
exp.frame.small <- data.frame(h = hs.small, mu = mus.small, cross.entropy = c(as.vector(results.array[,2,]), as.vector(results.array.para[,2]) ,as.vector(results.array.ml[,2]), rep(ce, length(h.set)) ) )

# model selection index, exp frame ???
# h = 1.1
# mu = 0.0008


exp.frame.small <- exp.frame.small[exp.frame.small$mu %in% c("0","0.00028","8e-04","0.01549","0.3","Parametric","Unregularized","Z-Score Matching"), ]


exp.frame.small$mu <- factor(exp.frame.small$mu, levels = c("0","0.00028","8e-04","0.01549","0.3","Parametric","Unregularized","Z-Score Matching"))
library(RColorBrewer)
library(colorRamps)
library(pals)

library(Polychrome)

# build-in color pal



P5 <- brewer.pal(6, "Blues")
P5 <- P5[2:6]
P8 <- c(P5,"green","red","orange")

P12 <- c(brewer.pal(10, "Blues")[3:10],"green","red","orange","purple")
gauss.plot <- ggplot(data = gauss.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + geom_line()  + ggtitle("Gaussian Kernel")  + ylab("Crosswalk Sample Cross Entropy") + coord_cartesian(xlim = c(0,4), ylim = c(950, 1100)) + scale_color_manual(values= as.vector(P12)) + annotate("point", x = gauss.frame$h[48], y = gauss.frame$cross.entropy[48], colour = "red", shape = "x", size = 4) + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                                                                                                                                                                                                       axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                                                                                                                                                                                                       axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                                                                                                                       axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                                                                                                                       title = element_text( size = title.size),
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.text = element_text(size=axis.ticks.size))

exp.plot <- ggplot(data = exp.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + geom_line()  + ggtitle("Laplace Kernel")  + ylab("Crosswalk Sample Cross Entropy")  + coord_cartesian(xlim = c(0,4), ylim = c(950, 1100))  + scale_color_manual(values= as.vector(P12)) + annotate("point", x = exp.frame$h[29], y = exp.frame$cross.entropy[29], colour = "blue", shape = "x", size = 4)+ theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                                                                                                                                                                                                axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                                                                                                                                                                                                axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                                                                                                                axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                                                                                                                title = element_text( size = title.size),
                                                                                                                                                                                                                                                                                                                                                                                                                legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                                                                                                                                                                                                legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                                                                                                                                                                                                legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                                                                                                                                                                                                legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                                                                                                                                                                                                legend.text = element_text(size=axis.ticks.size))

triangle.plot <- ggplot(data = triangle.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + geom_line()  + ggtitle("Triangle Kernel")  + ylab("Crosswalk Sample Cross Entropy")  + coord_cartesian(xlim = c(0,12), ylim = c(950, 1100)) + scale_color_manual(values= as.vector(P12)) + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                                                                                         axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                                                                                         axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                         axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                         title = element_text( size = title.size),
                                                                                                                                                                                                                                                                                                         legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                                                                                         legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                                                                                         legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                                                                                         legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                                                                                         legend.text = element_text(size=axis.ticks.size))

epanechnikov.plot <- ggplot(data = epanechnikov.frame, aes(x = h, y = cross.entropy, color = mu, group = mu)) + geom_line()  + ggtitle("Epan. Kernel")  + ylab("Crosswalk Sample Cross Entropy")  + coord_cartesian(xlim = c(0,12), ylim = c(950, 1100)) + scale_color_manual(values= as.vector(P12)) + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                                                                                                     axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                                                                                                     axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                     axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                     title = element_text( size = title.size),
                                                                                                                                                                                                                                                                                                                     legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                                                                                                     legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                                                                                                     legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                                                                                                     legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                                                                                                     legend.text = element_text(size=axis.ticks.size))
# argmin is the gaussian kernel with mu = 0.01







exp.plot.small <- ggplot(data = exp.frame.small, aes(x = h, y = cross.entropy, color = mu, group = mu)) + geom_line()  + ggtitle("Laplace Kernel")  + ylab("Crosswalk Sample Cross Entropy")  + coord_cartesian(xlim = c(0,4), ylim = c(950, 1100))  + scale_color_manual(values= as.vector(P8)) + annotate("point", x = exp.frame$h[20], y = exp.frame$cross.entropy[20], colour = "blue", shape = "x", size = 8)+ theme(axis.text.x = element_text( size = axis.ticks.size),
                                                                                                                                                                                                                                                                                                                                                                                                                axis.text.y = element_text( size = axis.ticks.size),  
                                                                                                                                                                                                                                                                                                                                                                                                                axis.title.x = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                                                                                                                axis.title.y = element_text( size = axis.size),
                                                                                                                                                                                                                                                                                                                                                                                                                title = element_text( size = title.size),
                                                                                                                                                                                                                                                                                                                                                                                                                legend.key.size = unit(1, 'cm'), #change legend key size
                                                                                                                                                                                                                                                                                                                                                                                                                legend.key.height = unit(1, 'cm'), #change legend key height
                                                                                                                                                                                                                                                                                                                                                                                                                legend.key.width = unit(1, 'cm'), #change legend key width
                                                                                                                                                                                                                                                                                                                                                                                                                legend.title = element_text(size=axis.size), #change legend title font size
                                                                                                                                                                                                                                                                                                                                                                                                                legend.text = element_text(size=axis.ticks.size))


exp.plot.small


 
gauss.plot  
exp.plot 
triangle.plot 
epanechnikov.plot

grid.arrange(arrangeGrob(gauss.plot, exp.plot, ncol = 2),
             arrangeGrob(triangle.plot, epanechnikov.plot, ncol = 2),
             nrow = 2)   



png(filename = "plots/conversion_NACC_MMSE_MOCA.png",
    width = png.width, height = png.height)

exp.plot.small

# Close the pdf file
dev.off() 


png(filename = "plots/conversion_NACC_MMSE_MOCA_large.png",
    width = png.width, height = png.height)


ggarrange(gauss.plot, exp.plot,
          triangle.plot, epanechnikov.plot,
          ncol=2, nrow=2, common.legend = TRUE, legend="right") 

# Close the pdf file
dev.off() 


min(epanechnikov.frame$cross.entropy, na.rm = T)










###### AAIC Small Plot 

cw.test <- read.csv(file = "data/NACCMMSE_to_MOCATOTS_test.csv")

# We now test the conversion error 
opt.mu.y <- 0.002
opt.mu.z <- 0.0008

mu.y <- opt.mu.y # 0.002
mu.z <- opt.mu.z # 0.0008


# Assume we have already computed an optimal h ker and mu from the first model 
Ny = 30
Nz = 30

num.bins <- 1000 #300

ker.z = exponential_kernel
hz = 1.1

A.model.y <- A.matrix.binomial(R_bins = num.bins, N = Ny, 
                               numeric.points = 400)

A.model.z <- A.matrix.compute(R_bins = num.bins, N = Nz, 
                              ker = exponential_kernel, h = 1.1,
                              numeric.points = 400)

h.set<- c(0.35,1.0,hz,1.5,2.0,2.5,4.0,6.0,8.0,12.0) 
results.array <- array(NA, dim = c(length(h.set), 4))


### population test set 


age.range <- 60:85
group.range <- 1:4
x.design <- data.frame(age = rep(age.range, times = length(group.range)), 
                       group = rep(group.range, each = length(age.range)))

p.hat.y.set <- cond.dist.est(x.design = x.design, train.data = y.dat, outcome = "y", N = Ny, age.window = 3)
p.hat.z.set <- cond.dist.est(x.design = x.design, train.data = z.dat, outcome = "z", N = Nz, age.window = 3)


wide.p.hat.y <- cbind(x.design,p.hat.y.set)
colnames(wide.p.hat.y) <- c("age", "group", 0:Ny)

wide.p.hat.z <- cbind(x.design,p.hat.z.set)
colnames(wide.p.hat.z) <- c("age", "group", 0:Nz)



wide.p.hat.y$model.idx = 1:nrow(wide.p.hat.y)
wide.p.hat.z$model.idx = 1:nrow(wide.p.hat.z)

#### matching rows of the validation set to the observed data conditional estimate 

cw.test <- cw.test %>% left_join(wide.p.hat.y[,c('age', 'group','model.idx')], by = c('age', 'group')) 


hist.cols.y <- as.character(0:Ny)
hist.cols.z <- as.character(0:Nz)

latent.mix.list.short.y <- list()
model.observed.list.short.y <- list()
latent.mix.list.short.y.ml <- list()
model.observed.list.short.y.ml <- list()


for(m in 1:nrow(wide.p.hat.y)){
  train.p.hat.y <- as.numeric(wide.p.hat.y[m,hist.cols.y])
  model.estimate.y <- numeric.latent.fit(p.hat = train.p.hat.y, 
                                         A.matrix = A.model.y, 
                                         mu = opt.mu.y, show.plot = F)
  
  latent.mix.list.short.y[[m]] <- model.estimate.y$latent
  model.observed.list.short.y[[m]] <- model.estimate.y$observed
  
  model.estimate.y.ml <- numeric.latent.fit(p.hat = train.p.hat.y, 
                                            A.matrix = A.model.y, 
                                            mu = 0, show.plot = F)
  
  latent.mix.list.short.y.ml[[m]] <- model.estimate.y.ml$latent
  model.observed.list.short.y.ml[[m]] <- model.estimate.y.ml$observed

  
  cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat.y)), end="\r")
}

cat(end="\n")

latent.mix.list.y <- list()
model.observed.list.y <- list()

latent.mix.list.y.ml <- list()
model.observed.list.y.ml <- list()

for(m in 1:nrow(cw.test)){
  k = cw.test$model.idx[m]
  latent.mix.list.y[[m]] <-latent.mix.list.short.y[[k]]
  model.observed.list.y[[m]] <- model.observed.list.short.y[[k]]
  
  latent.mix.list.y.ml[[m]] <- latent.mix.list.short.y.ml[[k]]
  model.observed.list.y.ml[[m]] <- model.observed.list.short.y.ml[[k]]
}



### Fixed Model Estimate of y varying z 

for(i in 1:length(h.set)){
  hz.tmp <- h.set[i]
  A.model.z <- A.matrix.compute(R_bins = num.bins, N = Nz, 
                                ker = exponential_kernel, h = hz.tmp,
                                numeric.points = 400)
  
  latent.mix.list.short.z <- list()
  model.observed.list.short.z <- list()
  latent.mix.list.short.z.ml <- list()
  model.observed.list.short.z.ml <- list()
  
  
  for(m in 1:nrow(wide.p.hat.z)){
    train.p.hat.z <- as.numeric(wide.p.hat.z[m,hist.cols.z])
    model.estimate.z <- numeric.latent.fit(p.hat = train.p.hat.z, 
                                           A.matrix = A.model.z, 
                                           mu = opt.mu.z, show.plot = F)
    
    latent.mix.list.short.z[[m]] <- model.estimate.z$latent
    model.observed.list.short.z[[m]] <- model.estimate.z$observed
    
    model.estimate.z.ml <- numeric.latent.fit(p.hat = train.p.hat.z, 
                                              A.matrix = A.model.z, 
                                              mu = 0, show.plot = F)
    
    latent.mix.list.short.z.ml[[m]] <- model.estimate.z.ml$latent
    model.observed.list.short.z.ml[[m]] <- model.estimate.z.ml$observed
    
    
    cat(paste0("Latent Model Computed Row: ", m,"/",nrow(wide.p.hat.z)), end="\r")
  }
  
  latent.mix.list.z <- list()
  model.observed.list.z <- list()
  
  latent.mix.list.z.ml <- list()
  model.observed.list.z.ml <- list()
  
  for(m in 1:nrow(cw.test)){
    k = cw.test$model.idx[m]
    latent.mix.list.z[[m]] <-latent.mix.list.short.z[[k]]
    model.observed.list.z[[m]] <- model.observed.list.short.z[[k]]
    
    latent.mix.list.z.ml[[m]] <- latent.mix.list.short.z.ml[[k]]
    model.observed.list.z.ml[[m]] <- model.observed.list.short.z.ml[[k]]
  }
  
  
  
  res <- score.conversion(test.pairs = cw.test[,c(1,2)], 
                          latent.mix.list.y = latent.mix.list.y,
                          model.observed.list.y = model.observed.list.y,
                          latent.mix.list.z = latent.mix.list.z, 
                          model.observed.list.z = model.observed.list.z, 
                          Ny = Ny, 
                          Nz = Nz, ker.z = ker.z, hz = hz.tmp, 
                          n.samp = 4, method = "CrossEntropy", binomial.y = T)
  
  res.ml <- score.conversion(test.pairs = cw.test[,c(1,2)], 
                             latent.mix.list.y = latent.mix.list.y.ml,
                             model.observed.list.y = model.observed.list.y.ml,
                             latent.mix.list.z = latent.mix.list.z.ml, 
                             model.observed.list.z = model.observed.list.z.ml, 
                             Ny = Ny, 
                             Nz = Nz, ker.z = ker.z, hz = hz.tmp, 
                             n.samp = 4, method = "CrossEntropy", binomial.y = T)
  
  results.array[i,3] <- res.ml
  results.array[i,4] <- res
  
}



muy <- mean(y.dat$y)
sdy <- sd(y.dat$y)
muz <- mean(z.dat$z)
sdz <- sd(z.dat$z)



ce.zscore = 0

for(i in 1:nrow(cw.test)){
  z.prob = naive.conversion.prob(y = cw.test[i,1], muy = muy, sdy = sdy,
                                 muz = muz, sdz = sdz, Nz = Nz)
  ce.zscore = ce.zscore - log(z.prob[cw.test[i,2] + 1])
  cat(paste0("Test sample ", i, "/",nrow(cw.test)), end = "\r")
}
cat(end = "\n")




# quantile matching method: 

# naive method: 

quantile_map_yz <- c()
for(y in 0:Ny){
  cdf.y <- mean(  y.dat$y <= y )
  z.pred <- round(quantile(z.dat$z, cdf.y))
  quantile_map_yz[y+1] <- z.pred
}







ce.quant.match = 0

for(i in 1:nrow(cw.test)){
  z.prob = quantile.conversion.prob(y = cw.test[i,1],quantile_map_yz = quantile_map_yz, decay_rate = 0.25, Nz = Nz)
  ce.quant.match = ce.quant.match -log(z.prob[cw.test[i,2] + 1])
  cat(paste0("Test sample ", i, "/",nrow(cw.test)), end = "\r")
}

results.array[,1] <- ce.zscore
results.array[,2] <- ce.quant.match



ce.vec <- c(results.array[,1],results.array[,2],results.array[,3],results.array[,4])

h.vec <- rep(h.set, 4)
method.vec <- rep(c("Z Score", "Quantile Matching", "Nonparametric", "Reg. Nonparametric"), each = length(h.set))
method.vec <- factor(method.vec, levels = c("Z Score", "Quantile Matching", "Nonparametric", "Reg. Nonparametric"))

results.data <- data.frame(h.vec = h.vec, ce.vec = ce.vec, method = method.vec)
#saveRDS(results.data,"data/AAIC_MMSE_MOCA_plot_data.rds")
results.data1 <- readRDS("data/AAIC_MMSE_MOCA_plot_data.rds")
hz = 1.1
aaic_results_plot1 <- ggplot(data = results.data1, aes(x = h.vec, y = ce.vec, group = method, color = method)) + 
  geom_line() + ylab("Conversion Cross Entropy") + xlab("Bandwidth (h)") + ggtitle("MMSE to MOCA") + 
  annotate("point", x = hz, y = results.data1$ce.vec[results.data1$h.vec == hz & results.data1$method == "Reg. Nonparametric"], 
           colour = "blue", shape = "x", size = 8) + 
  coord_cartesian(xlim = c(0,12), ylim = c(900, 1100)) + theme(axis.text.x = element_text( size = axis.ticks.size),
                                                               axis.text.y = element_text( size = axis.ticks.size),  
                                                               axis.title.x = element_text( size = axis.size),
                                                               axis.title.y = element_text( size = axis.size),
                                                               title = element_text( size = title.size),
                                                               legend.key.size = unit(1, 'cm'), #change legend key size
                                                               legend.key.height = unit(1, 'cm'), #change legend key height
                                                               legend.key.width = unit(1, 'cm'), #change legend key width
                                                               legend.title = element_text(size=axis.size), #change legend title font size
                                                               legend.text = element_text(size=axis.ticks.size)) + 
  labs(col="Method")

aaic_results_plot1




















#### Parametric Approximation: 

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
############
latent.variational.moments <- renormalized_r_function(N, is.binomial = T)
library(caret)

y.dat$group <- factor(y.dat$group, levels = c(1,2,3,4))
dmy <- dummyVars(" ~ .", data = y.dat)
y.dat.wide <- data.frame(predict(dmy, newdata = y.dat))
y.dat.wide


z.dat$group <- factor(z.dat$group, levels = c(1,2,3,4))
dmy <- dummyVars(" ~ .", data = z.dat)
z.dat.wide <- data.frame(predict(dmy, newdata = z.dat))
z.dat.wide


params.y <- fit_logitnorm_model(y.dat.wide$y, X = y.dat.wide[,2:6], latent.variational.moments)
params.z <- fit_logitnorm_model(z.dat.wide$z, X = z.dat.wide[,2:6], latent.variational.moments)


cw.test$group <- factor(cw.test$group, levels = c(1,2,3,4))
dmy <- dummyVars(" ~ .", data = cw.test)
cw.test.wide <- data.frame(predict(dmy, newdata = cw.test))

############

logistic <- function(x){
  return(1/(1 + exp(-x)))
}

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



# parametric sample as in the non-parametric version. I.e. sample normal and logit transform 

y = y
z = z
x = x.row
params.y = params.y
params.z = params.z
Ny = Ny
Nz = Nz
ker.z = ker.tmp
hz = h.tmp
binomial.y = T
R_bins = 1000

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




## Estimated parameters 
params.y
params.z





















