
# analysis

load("results.Rdata")

sim_number <- length(reg_result)

##############################################
# compare the p values

# individual level
p_individual <- array(dim=c(sim_number,1+12))

for(i in 1:sim_number){
  p_individual[i,] <- reg_result[[i]]$individual$p
}

apply(ifelse(p_individual<0.05,1,0),2,sum)/sim_number

# the factor levels
N <- rep(0,sim_number)
for(i in 1:sim_number){
  N[i] <- reg_result[[i]]$N
}

table(N)

# look at 1-3 PC factors

# 1-factor

p_1 <- array(dim=c(sim_number,1+1))

for(i in 1:sim_number){
  p_1[i,] <- reg_result[[i]]$principal[[1]]$p
}

apply(ifelse(p_1<0.05,1,0),2,sum)/sim_number

# 2-factor

# check for abnormal loadings
loadings <- array(dim=c(sim_number,12,2))
for(i in 1:sim_number){
  loadings[i,,] <- reg_result[[i]]$principal[[2]]$loadings
}

loadings <- ifelse(abs(loadings)>0.5,1,0)

apply(loadings,c(2,3),sum)
# factor 1 - (1,2); factor 2 - (3)

p_2 <- array(dim=c(sim_number,1+2))

for(i in 1:sim_number){
  p_2[i,] <- reg_result[[i]]$principal[[2]]$p
}

apply(ifelse(p_2<0.05,1,0),2,sum)/sim_number

# 3-factor

# check for abnormal loadings
loadings <- array(dim=c(sim_number,12,3))
for(i in 1:sim_number){
  if(length( reg_result[[i]]$principal)>2){
    loadings[i,,] <- reg_result[[i]]$principal[[3]]$loadings
  }
}

loadings <- ifelse(abs(loadings)>0.5,1,0)

# the factors change by every simulation

apply(loadings,c(2,3),function(x){sum(x,na.rm=T)})

# do it this way: multiply by loadings
p_3 <- array(dim=c(sim_number,12))

for(i in 1:sim_number){
  if(length( reg_result[[i]]$principal)>2){
    for(j in 1:12){
      p_3[i,j] <- sum(reg_result[[i]]$principal[[3]]$p[2:4]*loadings[i,j,])
    }
  }  
}

apply(ifelse(p_3<0.05,1,0),2,function(x){mean(x,na.rm=T)})


















