
library(stats)
library(psych)
library(dendextend)
library(foreign)
library(lavaan)
library(semTools)
library(MBESS)
library(GPArotation)

require(MuMIn)
require(Metrics)
require(rstanarm)
require(loo)
require(bridgesampling)

require(psych)

# 500

set.seed(117)

load("parameters.Rdata")

res <- list()

sim_number <- 1000

for(l in 1:sim_number){
  
  skip_to_next <- FALSE

  k <- sample(1:100,1)
  
  fx <- parameters$fx_v2[[sample(1:10000,1)]]
  fy <- parameters$fy
  Phi <- parameters$Phi_v[[k]]
  
  sample_size <- 500
  
  aic <- rep(0,4)
  bic <- rep(0,4)
  bf <- list()
  
  # data
  bifdat <- sim.structure(fx = fx, Phi = Phi, fy=fy, n = sample_size)$observed
  dat <- as.data.frame(bifdat)
  
  # individual
  fit1 <- lm(Vy1 ~ Vx1 + Vx2 + Vx3 + Vx4 + Vx5 + Vx6 + Vx7 + Vx8 + Vx9 + Vx10 + Vx11 + Vx12, data=dat)
  aic[1] <- AIC(fit1)
  bic[1] <- BIC(fit1)
  
  warns.b <- tryCatch(fit1.b <- stan_glm(Vy1 ~ Vx1 + Vx2 + Vx3 + Vx4 + Vx5 + Vx6 + Vx7 + Vx8 + Vx9 + Vx10 + Vx11 + Vx12,
                                         data = dat, prior = normal(location=rep(0,12),scale=rep(1,12)),chains=1,
                                         diagnostic_file = file.path(tempdir(), "df.csv")),error=function(e) { skip_to_next <<- TRUE}, warning=function(w) w)
  warns.bf <- tryCatch(bf[[1]] <- bridge_sampler(fit1.b),error=function(e) { skip_to_next <<- TRUE}, warning=function(w) w)
  
  if(skip_to_next) { next } 
  
  # psych
  warns <- tryCatch(fa.parallel(dat[,1:12]),error=function(e) { skip_to_next <<- TRUE}, warning=function(w) w)
  
  # level-1
  fa <- pca(dat[,1:12],nfactors=1)
  scores <- fa$scores
  fit2 <- lm(Vy1 ~ scores, data=dat)
  aic[2] <- AIC(fit2)
  bic[2] <- BIC(fit2)
  
  fit2.b <- stan_glm(Vy1 ~ scores,
                     data = dat, prior = normal(location=rep(0,1),scale=rep(1,1)),chains=1,
                     diagnostic_file = file.path(tempdir(), "df.csv"))
  bf[[2]] <- bridge_sampler(fit2.b)
  
  # level-2
  fa <- pca(dat[,1:12],nfactors=2)
  scores <- fa$scores
  weights <- c(1,1,1,1,1,1,1,1,2,2,2,2)
  w2 <- c(sum(abs(fa$loadings[,1])*weights)/4,
          sum(abs(fa$loadings[,2])*weights)/4)
  w2 <- rank(w2)
  fit3 <- lm(Vy1 ~ scores, data=dat)
  aic[3] <- AIC(fit3)
  bic[3] <- BIC(fit3)
  
  fit3.b <- stan_glm(Vy1 ~ scores,
                     data = dat, prior = normal(location=rep(0,2),scale=rep(1,2)),chains=1,
                     diagnostic_file = file.path(tempdir(), "df.csv"))
  bf[[3]] <- bridge_sampler(fit3.b)
  
  # level-3
  fa <- pca(dat[,1:12],nfactors=3)
  scores <- fa$scores
  weights <- c(1,1,1,1,2,2,2,2,3,3,3,3)
  w3 <- c(sum(abs(fa$loadings[,1])*weights)/4,
          sum(abs(fa$loadings[,2])*weights)/4,
          sum(abs(fa$loadings[,3])*weights)/4)
  w3 <- rank(w3)
  fit4 <- lm(Vy1 ~ scores, data=dat)
  aic[4] <- AIC(fit4)
  bic[4] <- BIC(fit4)
  
  fit4.b <- stan_glm(Vy1 ~ scores,
                     data = dat, prior = normal(location=rep(0,3),scale=rep(1,3)),chains=1,
                     diagnostic_file = file.path(tempdir(), "df.csv"))
  bf[[4]] <- bridge_sampler(fit4.b)
  
  
  fit <- list(fit1,fit2,fit3,fit4,fit1.b,fit2.b,fit3.b,fit4.b)
  
  res[[l]] <- list(aic=match(min(aic),aic),
                   p_aic=round(summary(fit[[match(min(aic),aic)]])$coefficients[,4],digits=3),
                   bic=match(min(bic),bic),
                   p_bic=round(summary(fit[[match(min(bic),bic)]])$coefficients[,4],digits=3),
                   k=k,
                   Phi=Phi,
                   fx=fx,
                   warning=warns$message,
                   nfact=warns$nfact,
                   w2=w2,
                   w3=w3,
                   bf=bf)
  
  print(l)
}


save(res,file="model_results_condition_7_500.Rdata")

