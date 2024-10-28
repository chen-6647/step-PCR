
library(stats)
library(psych)
library(dendextend)
library(foreign)
library(lavaan)
library(semTools)
library(MBESS)
library(GPArotation)


fx <- matrix(c(0.4, 0.6, 0, 0,
               0.4, 0.6, 0, 0,
               0.4, 0.6, 0, 0,
               0.4, 0.6, 0, 0,
               0.4, 0, 0.6, 0,
               0.4, 0, 0.6, 0,
               0.4, 0, 0.6, 0,
               0.4, 0, 0.6, 0,
               0.4, 0, 0, 0.6,
               0.4, 0, 0, 0.6,
               0.4, 0, 0, 0.6,
               0.4, 0, 0, 0.6),12,4,byrow = TRUE) 
fy <- matrix(1,1,1,byrow = TRUE) 
Phi <- matrix(c(1, 0, 0, 0, 0.6,
                0, 1, 0.6, 0, 0,
                0, 0.6, 1, 0, 0,
                0, 0, 0, 1, 0,
                0.6, 0, 0, 0, 1), 5,5, byrow = TRUE) 



###############################################################
# AIC, BIC, R2, R2_adj, LOOIC, WAIC, BF
###############################################################

# sample size: 100, 800
# sim_number: 200

sample_size <- 200
sim_number <- 1000

# three parts
# (1) parallel analysis
# (2) principal component and regression outcomes
# (3) model comparison and selection indices

# frequentist methods
# AIC: AIC(fit)
# BIC: BIC(fit)
# R2: summary(fit)$r.squared
# R2_adj: summary(fit)$adj.r.squared

# Bayesian methods
# use rstanarm and loo
# looic: loo(fit.b) -> a$looic (will be depreciated)
# waic: waic(fit.b) -> a$waic (will be depreciated)
# Bayes factor: package "bridgesampling" 


require(MuMIn)
require(Metrics)
require(rstanarm)
require(loo)
require(bridgesampling)

data <- list()

set.seed(117)

for(k in 1:sim_number){
  
  bifdat <- sim.structure(fx = fx, Phi = Phi, fy=fy, n = sample_size)$observed
  bifdat <- as.data.frame(bifdat)
  
  data[[k]] <- bifdat
  
}

save(data,file="data.Rdata")


########################################
########################################

# load saved data and calculate the principal component, regression, and model comparison outcomes


load("data.Rdata")

sim_number <- length(data)
sample_size <- dim(data[[1]])[1]

reg_result <- list()

for(i in 1:sim_number){
  dat <- data[[i]]
  pal <- fa.parallel(dat[,1:12])
  res <- list(ncomp=pal$ncomp,nfact=pal$nfact)
  
  N <- res$nfact  # the maximum number of PCs to investigate
  
  # regressions
  # model 1 - individual level
  # frequentist methods
  fit <- lm(Vy1 ~ Vx1 + Vx2 + Vx3 + Vx4 + Vx5 + Vx6 + Vx7 + Vx8 + Vx9 + Vx10 + Vx11 + Vx12, data=dat)
  p <- summary(fit)$coefficients[,4]
  aic <- AIC(fit)
  bic <- BIC(fit)
  r2 <- summary(fit)$r.squared
  r2_adj <- summary(fit)$adj.r.squared
  
  # Bayesian methods - use N(0,1) priors
  fit.b <- stan_glm(Vy1 ~ Vx1 + Vx2 + Vx3 + Vx4 + Vx5 + Vx6 + Vx7 + Vx8 + Vx9 + Vx10 + Vx11 + Vx12,
                   data = dat, prior = normal(location=rep(0,12),scale=rep(1,12)),chains=1,
                   diagnostic_file = file.path(tempdir(), "df.csv"))
  
  looic <- suppressWarnings(loo(fit.b)$looic)
  waic <- suppressWarnings(waic(fit.b)$waic)
  bf <- bridge_sampler(fit.b)
  
  individual <- list(p=p,
                     aic=aic,
                     bic=bic,
                     r2=r2,
                     r2_adj=r2_adj,
                     looic=looic,
                     waic=waic,
                     bf=bf)
  
  
  # PCR models
  principal <- list()
  
  for(n in 1:N){
    fa <- pca(dat[,1:12],nfactors=n)
    scores <- fa$scores
    
    # frequentist
    fit <- lm(Vy1 ~ scores, data=dat)
    p <- summary(fit)$coefficients[,4]
    aic <- AIC(fit)
    bic <- BIC(fit)
    r2 <- summary(fit)$r.squared
    r2_adj <- summary(fit)$adj.r.squared
    
    # Bayesian methods - use N(0,1) priors
    fit.b <- stan_glm(Vy1 ~ scores,
                      data = dat, prior = normal(location=rep(0,n),scale=rep(1,n)),chains=1,
                      diagnostic_file = file.path(tempdir(), "df.csv"))
    
    looic <- suppressWarnings(loo(fit.b)$looic)
    waic <- suppressWarnings(waic(fit.b)$waic)
    bf <- bridge_sampler(fit.b)
    
    
    principal[[n]] <- list(N_pc=n,
                           loadings=fa$loadings,
                           p=p,
                           aic=aic,
                           bic=bic,
                           r2=r2,
                           r2_adj=r2_adj,
                           looic=looic,
                           waic=waic,
                           bf=bf)
  }
  
  reg_result[[i]] <- list(ncomp=pal$ncomp,
                          nfact=pal$nfact,
                          N=N,
                          individual=individual,
                          principal=principal)
}




save(reg_result,file="results.Rdata")



