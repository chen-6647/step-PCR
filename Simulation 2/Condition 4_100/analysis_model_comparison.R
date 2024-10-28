

require(rstanarm)
require(loo)
require(bridgesampling)

# analysis

load("results.Rdata")

sim_number <- length(reg_result)

# summarize the model comparison results

N <- rep(0,sim_number)
for(i in 1:sim_number){
  N[i] <- reg_result[[i]]$N
}

N <- max(N)

get_best <- function(reg_result,sim_number,N,var.name,min=T){
  
  index <- array(NA,dim=c(sim_number,1+N))
  
  for(i in 1:sim_number){
    
    index[i,1] <- reg_result[[i]]$individual[[var.name]]
    
    for(n in 1:reg_result[[i]]$N){
      res <- reg_result[[i]]$principal[[n]]
      index[i,1+n] <- res[[var.name]]
    }
  }
  if(min==T){
    best <- apply(index,1,function(x){match(min(x,na.rm=T),x)})   
  }
  else{
    best <- apply(index,1,function(x){match(max(x,na.rm=T),x)})    
  }
  return(best)
}

find_best_bf <- function(reg_result,sim_number){
  best <- rep(0,sim_number)
  
  for(i in 1:sim_number){
    bf <- list()
    bf[[1]] <- reg_result[[i]]$individual$bf
    for(n in 1:reg_result[[i]]$N){
      bf[[n+1]] <- reg_result[[i]]$principal[[n]]$bf
    }
    k <- ifelse(bf(bf[[1]],bf[[2]])[1]>1,1,2)
    for(n in 2:(reg_result[[i]]$N)){
      k <- ifelse(bf(bf[[k]],bf[[n+1]])[1]>1,k,n+1)
    }
    best[i] <- k
  }
  return(best)
}


mc.results <- data.frame(array(dim=c(N+1,7)))

names(mc.results) <- c("aic","bic","r2","r2_adj","looic","waic","bf")
row.names(mc.results) <- 1:(N+1)

# aic
best <- get_best(reg_result,sim_number,N,var.name="aic",min=T)
for(n in 1:(N+1)){
  mc.results[n,"aic"] <- sum(ifelse(best==n,1,0))
}

# bic
best <- get_best(reg_result,sim_number,N,var.name="bic",min=T)
for(n in 1:(N+1)){
  mc.results[n,"bic"] <- sum(ifelse(best==n,1,0))
}

# r2
best <- get_best(reg_result,sim_number,N,var.name="r2",min=F)
for(n in 1:(N+1)){
  mc.results[n,"r2"] <- sum(ifelse(best==n,1,0))
}

# r2_adj
best <- get_best(reg_result,sim_number,N,var.name="r2_adj",min=F)
for(n in 1:(N+1)){
  mc.results[n,"r2_adj"] <- sum(ifelse(best==n,1,0))
}

# looic
best <- get_best(reg_result,sim_number,N,var.name="looic",min=T)
for(n in 1:(N+1)){
  mc.results[n,"looic"] <- sum(ifelse(best==n,1,0))
}

# waic
best <- get_best(reg_result,sim_number,N,var.name="waic",min=T)
for(n in 1:(N+1)){
  mc.results[n,"waic"] <- sum(ifelse(best==n,1,0))
}

# bf 
best <- find_best_bf(reg_result,sim_number)
for(n in 1:(N+1)){
  mc.results[n,"bf"] <- sum(ifelse(best==n,1,0))
}

write.csv(mc.results,file="model_comparison_results.csv")

