
# make a list of matrices: fx, fy, Phi

fy <- matrix(1,1,1,byrow = TRUE) 

# make 1000 types of fx

set.seed(117)

number <- 10000

fx <- list()

for(i in 1:number){
  x1 <- runif(12,0.3,0.5)
  x2 <- c(runif(4,0.5,0.8),runif(4,0,0.3),runif(4,0,0.3))
  x3 <- c(runif(4,0,0.3),runif(4,0.5,0.8),runif(4,0,0.3))
  x4 <- c(runif(4,0,0.3),runif(4,0,0.3),runif(4,0.5,0.8))
  fx[[i]] <- matrix(cbind(x1,x2,x3,x4),12,4,byrow = F)
}

Phi0 <- matrix(c(1, 0, 0, 0, 0,
                 0, 1, 0, 0, 0,
                 0, 0, 1, 0, 0,
                 0, 0, 0, 1, 0,
                 0, 0, 0, 0, 1), 5,5, byrow = TRUE) 

# 10 types each

Phi <- list()

# no correlations

for(i in 1:100){
  p <- Phi0
  p[2,3] <- p[3,2] <- runif(1,0.3,0.6)
  Phi[[i]] <- p
}

for(i in 101:200){
  p <- Phi0
  p[2,3] <- p[3,2] <- runif(1,0.3,0.6)
  p[1,5] <- p[5,1] <- runif(1,0.3,0.6)
  Phi[[i]] <- p  
}

for(i in 201:300){
  p <- Phi0
  p[2,3] <- p[3,2] <- runif(1,0.3,0.6)
  p[2,5] <- p[5,2] <- runif(1,0.2,0.4)
  p[3,5] <- p[5,3] <- runif(1,0.2,0.4)
  Phi[[i]] <- p   
}

for(i in 301:400){
  p <- Phi0
  p[2,3] <- p[3,2] <- runif(1,0.3,0.6)
  p[2,5] <- p[5,2] <- runif(1,0.3,0.6)
  Phi[[i]] <- p   
}

for(i in 401:500){
  p <- Phi0
  p[2,3] <- p[3,2] <- runif(1,0.3,0.6)
  p[4,5] <- p[5,4] <- runif(1,0.3,0.6)
  Phi[[i]] <- p   
}



fx_v1 <- list()

for(i in 1:number){
  x1 <- runif(12,0.3,0.5)
  x2 <- c(runif(4,0.5,0.8),runif(4,0,0.3),runif(4,0,0.3))
  x3 <- c(runif(4,0,0.3),runif(4,0.5,0.8),runif(4,0,0.3))
  x4 <- c(runif(4,0,0.3),runif(4,0,0.3),runif(4,0.5,0.8))
  x5 <- c(runif(1,0.2,0.6),rep(0,11))
  fx_v1[[i]] <- matrix(cbind(x1,x2,x3,x4,x5),12,5,byrow = F)
}

fx_v2 <- list()

for(i in 1:number){
  x1 <- runif(12,0.3,0.5)
  x2 <- c(runif(4,0.5,0.8),runif(4,0,0.3),runif(4,0,0.3))
  x3 <- c(runif(4,0,0.3),runif(4,0.5,0.8),runif(4,0,0.3))
  x4 <- c(runif(4,0,0.3),runif(4,0,0.3),runif(4,0.5,0.8))
  x5 <- c(rep(0,8),runif(1,0.2,0.6),rep(0,3))
  fx_v2[[i]] <- matrix(cbind(x1,x2,x3,x4,x5),12,5,byrow = F)
}

Phi0_v <- matrix(c(1, 0, 0, 0, 0, 0,
                   0, 1, 0, 0, 0, 0,
                   0, 0, 1, 0, 0, 0,
                   0, 0, 0, 1, 0, 0,
                   0, 0, 0, 0, 1, 1,
                   0, 0, 0, 0, 1, 1), 6, 6, byrow = TRUE) 

Phi_v <- list()

# no correlations

for(i in 1:100){
  p <- Phi0_v
  p[2,3] <- p[3,2] <- runif(1,0.3,0.6)
  Phi_v[[i]] <- p
}


parameters <- list(fx=fx,fy=fy,Phi=Phi,fx_v1=fx_v1,fx_v2=fx_v2,Phi_v=Phi_v)

save(parameters,file="parameters.Rdata")




