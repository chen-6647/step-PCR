
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

# Condition 1: no correlations

load("model_results_condition_7_500.Rdata")

res <- res[-which(sapply(res,is.null))]

sim_number <- length(res)

select <- c(4,4)

plot_name <- "Condition 7: Y ~ V9"


###########
# list how many times each level was selected
# and how many tests are statistically significant

aic_result <- list(m1=0,
                   p1=rep(0,12),
                   m2=0,
                   p2=rep(0,1),
                   m3=0,
                   p3=rep(0,2),
                   m4=0,
                   p4=rep(0,3))

bic_result <- list(m1=0,
                   p1=rep(0,12),
                   m2=0,
                   p2=rep(0,1),
                   m3=0,
                   p3=rep(0,2),
                   m4=0,
                   p4=rep(0,3))

bf_result <- list(m1=0,
                  p1=rep(0,12),
                  m2=0,
                  p2=rep(0,1),
                  m3=0,
                  p3=rep(0,2),
                  m4=0,
                  p4=rep(0,3))



for(l in 1:sim_number){
  if(is.null(res[[l]]$warning)){
    
    # AIC
    aic_result[[paste0("m",res[[l]]$aic)]] <- aic_result[[paste0("m",res[[l]]$aic)]] + 1
    if(res[[l]]$aic==1){
      if(length(res[[l]]$p_aic)==13){
        aic_result$p1 <- aic_result$p1 + ifelse(res[[l]]$p_aic[2:13]<.05,1,0)          
      }
    }
    if(res[[l]]$aic==2){
      aic_result$p2 <- aic_result$p2 + ifelse(res[[l]]$p_aic[2]<.05,1,0)
    }
    if(res[[l]]$aic==3){
      aic_result$p3 <- aic_result$p3 + ifelse(res[[l]]$p_aic[2:3]<.05,1,0)[match(c(1,2),res[[l]]$w2)]
    }
    if(res[[l]]$aic==4){
      aic_result$p4 <- aic_result$p4 + ifelse(res[[l]]$p_aic[2:4]<.05,1,0)[match(c(1,2,3),res[[l]]$w3)]
    }
    
    # BIC
    bic_result[[paste0("m",res[[l]]$bic)]] <- bic_result[[paste0("m",res[[l]]$bic)]] + 1
    if(res[[l]]$bic==1){
      if(length(res[[l]]$p_bic)==13){
        bic_result$p1 <- bic_result$p1 + ifelse(res[[l]]$p_bic[2:13]<.05,1,0)
      }
    }
    if(res[[l]]$bic==2){
      bic_result$p2 <- bic_result$p2 + ifelse(res[[l]]$p_bic[2]<.05,1,0)
    }
    if(res[[l]]$bic==3){
      bic_result$p3 <- bic_result$p3 + ifelse(res[[l]]$p_bic[2:3]<.05,1,0)[match(c(1,2),res[[l]]$w2)]
    }
    if(res[[l]]$bic==4){
      bic_result$p4 <- bic_result$p4 + ifelse(res[[l]]$p_bic[2:4]<.05,1,0)[match(c(1,2,3),res[[l]]$w3)]
    }
    
    # Bayes Factor
    bf <- res[[l]]$bf
    bf_optimal <- ifelse(is.null(bf[[1]]),2,ifelse(bf(bf[[1]],bf[[2]])[1]>1,1,2))
    bf_optimal <- ifelse(bf(bf[[bf_optimal]],bf[[3]])[1]>1,bf_optimal,3)
    bf_optimal <- ifelse(bf(bf[[bf_optimal]],bf[[4]])[1]>1,bf_optimal,4)
    
    bf_result[[paste0("m",bf_optimal)]] <- bf_result[[paste0("m",bf_optimal)]] + 1
  }
}

aic_result

bic_result

# Plot 1: Model selection results
height <-  c(aic_result$m2,
             bic_result$m2,
             bf_result$m2,
             aic_result$m3,
             bic_result$m3,
             bf_result$m3,
             aic_result$m4,
             bic_result$m4,
             bf_result$m4,
             aic_result$m1,
             bic_result$m1,
             bf_result$m1)

height <- height*3/sum(height)

cols <- rep(c("lightblue","blue","darkblue"),4)

density <- rep(c(NA,NA,20),4)

b <- barplot(height=height,col=cols,
             xlab="Model",ylab="Frequence of Selection",
             main=plot_name,ylim=c(0,1.1),xlim=c(0,18),
             density=density)

axis(1,at=c(1.9,5.5,9.1,12.7),labels=c('1-comp PCR','2-comp PCR','3-comp PCR','12-var MR'))

text(x=b,y=height+0.03,labels=paste0(round(height,2)*100,"%"))

rect(xleft=b[3*select[1]-1,1]-1.8,xright=b[3*select[2]-1,1]+1.8,ybottom=0,ytop=1.1,border='red')

legend("topright", title="Model Selection",
       c("AIC","BIC","BF"), fill=c("lightblue","blue","darkblue"), density=c(NA,NA,20), cex=0.8)

# aic bic table

output <- array(dim=c(4,1))

sum <- aic_result$m1 + aic_result$m2 + aic_result$m3 + aic_result$m4

output[1,1] <- paste(round(aic_result$m2/sum,digits=4)*100,"",
                     "&",                   
                     paste0("G: ", round(aic_result$p2/aic_result$m2,digits=4)*100,""),
                     "&",
                     round(bic_result$m2/sum,digits=4)*100,"",
                     "&", 
                     paste0("G: ", round(bic_result$p2/bic_result$m2,digits=4)*100,""))

output[2,1] <- paste(round(aic_result$m3/sum,digits=4)*100,"",
                     "&",  
                     paste0("F1F2: ", round(aic_result$p3[1]/aic_result$m3,digits=4)*100,";"),
                     paste0("F3: ", round(aic_result$p3[2]/aic_result$m3,digits=4)*100,""),
                     "&",
                     round(bic_result$m3/sum,digits=4)*100,"",
                     "&",  
                     paste0("F1F2: ", round(bic_result$p3[1]/bic_result$m3,digits=4)*100,";"),
                     paste0("F3: ", round(bic_result$p3[2]/bic_result$m3,digits=4)*100,""))

output[3,1] <- paste(round(aic_result$m4/sum,digits=4)*100,"",
                     "&", 
                     paste0("F1: ", round(aic_result$p4[1]/aic_result$m4,digits=4)*100,";"),
                     paste0("F2: ", round(aic_result$p4[2]/aic_result$m4,digits=4)*100,";"),
                     paste0("F3: ", round(aic_result$p4[3]/aic_result$m4,digits=4)*100,""),
                     "&",
                     round(bic_result$m4/sum,digits=4)*100,"",
                     "&", 
                     paste0("F1: ", round(bic_result$p4[1]/bic_result$m4,digits=4)*100,";"),
                     paste0("F2: ", round(bic_result$p4[2]/bic_result$m4,digits=4)*100,";"),
                     paste0("F3: ", round(bic_result$p4[3]/bic_result$m4,digits=4)*100,""))

output[4,1] <- paste(round(aic_result$m1/sum,digits=4)*100,"",
                     "&", 
                     paste0("V1-V12: ", round(min(aic_result$p1/aic_result$m1),digits=4)*100,"-",
                            round(max(aic_result$p1/aic_result$m1),digits=4)*100,""),
                     "&",
                     round(bic_result$m1/sum,digits=4)*100,"",
                     "&", 
                     paste0("V1-V12: ", round(min(bic_result$p1/bic_result$m1),digits=4)*100,"-",
                            round(max(bic_result$p1/bic_result$m1),digits=4)*100,""))

write.csv(output,"Table_output_7.csv")
