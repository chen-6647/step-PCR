

require(psych)

data <- read.csv("data.csv",sep="\t")


####### the dark triad

dd <- data[,11:22]

# 4 factors, 2 components
fa.parallel(dd)

a <- pca(dd,nfactors=3)
a$loadings # clearly some cross loadings

######## the dark triad associated with items in the Hypersensitive Narcissism scale

sub <- data[,11:22]

sub$HSN <- apply(data[,1:10],1,sum)

cor(sub)
# correlations mostly in the 0.3-0.5 area

######## direct assessment

set.seed(117)

sample_number <- 1000

size <- 500

p1 <- array(NA,dim=c(sample_number,12))

p_a <- array(NA,dim=c(sample_number,3))

p_b <- array(NA,dim=c(sample_number,1))

p2 <- array(NA,dim=c(sample_number,1))

p3 <- array(NA,dim=c(sample_number,2))

w2 <- array(NA,dim=c(sample_number,2))

p4 <- array(NA,dim=c(sample_number,3)) 

w3 <- array(NA,dim=c(sample_number,3))

BIC <- array(dim=c(sample_number,4))


for(l in 1:sample_number){
  subsub <- sub[sample(1:53981,size,replace=T),]
  parallel <- fa.parallel(subsub[,1:12])
  
  # individual model
  fit1 <- lm(HSN ~ DDP1 + DDP2 + DDP3 + DDP4 + DDN1 + DDN2 + DDN3 + DDN4 + DDM1 + DDM2 + DDM3 + DDM4,
             data=subsub)
  p1[l,] <- summary(fit1)$coefficients[2:13,4]
  BIC[l,1] <- BIC(fit1)
  
  # the partialing model

  subsub$DDP <- apply(subsub[,paste0("DDP",1:4)],1,sum)
  subsub$DDN <- apply(subsub[,paste0("DDN",1:4)],1,sum)
  subsub$DDM <- apply(subsub[,paste0("DDM",1:4)],1,sum)
  fit_a <- lm(HSN ~ DDP + DDN + DDM, data=subsub)
  p_a[l,] <- summary(fit_a)$coefficients[2:4,4]
  
  # model with only M as the predictor
  fit_b <- lm(HSN ~ DDM, data=subsub)
  p_b[l,] <- summary(fit_b)$coefficients[2,4]
  
  # level 1
  fa <- pca(subsub[,1:12],nfactors=1)
  scores <- fa$scores
  fit2 <- lm(HSN ~ scores, data=subsub)
  p2[l,] <- summary(fit2)$coefficients[2,4]
  BIC[l,2] <- BIC(fit2)
  
  # level 2
  fa <- pca(subsub[,1:12],nfactors=2)
  w2[l,1] <- ifelse(mean(fa$loadings[1:4,1])>0.5 & mean(fa$loadings[9:12,1])>0.5 & mean(fa$loadings[5:8,1])<0.5,1,
                    ifelse(mean(fa$loadings[1:4,1])<0.5 & mean(fa$loadings[9:12,1])<0.5 & mean(fa$loadings[5:8,1])>0.5,2,NA))
  w2[l,2] <- ifelse(mean(fa$loadings[1:4,2])>0.5 & mean(fa$loadings[9:12,2])>0.5 & mean(fa$loadings[5:8,2])<0.5,1,
                    ifelse(mean(fa$loadings[1:4,2])<0.5 & mean(fa$loadings[9:12,2])<0.5 & mean(fa$loadings[5:8,2])>0.5,2,NA))
  scores <- fa$scores
  fit3 <- lm(HSN ~ scores, data=subsub)
  p3[l,] <- summary(fit3)$coefficients[w2[l,]+1,4]
  BIC[l,3] <- BIC(fit3)
  
  # level 3
  if(parallel$nfact<3){
    next
  }
  else{
     fa <- pca(subsub[,1:12],nfactors=3)
     w3[l,1] <- ifelse(mean(fa$loadings[1:4,1])>0.5 & mean(fa$loadings[9:12,1])<0.5 & mean(fa$loadings[5:8,1])<0.5,1,
                    ifelse(mean(fa$loadings[1:4,1])<0.5 & mean(fa$loadings[9:12,1])>0.5 & mean(fa$loadings[5:8,1])<0.5,3,
                           ifelse(mean(fa$loadings[1:4,1])<0.5 & mean(fa$loadings[9:12,1])<0.5 & mean(fa$loadings[5:8,1])>0.5,2,NA)))
     w3[l,2] <- ifelse(mean(fa$loadings[1:4,2])>0.5 & mean(fa$loadings[9:12,2])<0.5 & mean(fa$loadings[5:8,2])<0.5,1,
                    ifelse(mean(fa$loadings[1:4,2])<0.5 & mean(fa$loadings[9:12,2])>0.5 & mean(fa$loadings[5:8,2])<0.5,3,
                           ifelse(mean(fa$loadings[1:4,2])<0.5 & mean(fa$loadings[9:12,2])<0.5 & mean(fa$loadings[5:8,2])>0.5,2,NA)))
     w3[l,3] <- ifelse(mean(fa$loadings[1:4,3])>0.5 & mean(fa$loadings[9:12,3])<0.5 & mean(fa$loadings[5:8,3])<0.5,1,
                    ifelse(mean(fa$loadings[1:4,3])<0.5 & mean(fa$loadings[9:12,3])>0.5 & mean(fa$loadings[5:8,3])<0.5,3,
                           ifelse(mean(fa$loadings[1:4,3])<0.5 & mean(fa$loadings[9:12,3])<0.5 & mean(fa$loadings[5:8,3])>0.5,2,NA)))
     scores <- fa$scores
     fit4 <- lm(HSN ~ scores, data=subsub)
     p4[l,] <- summary(fit4)$coefficients[w3[l,]+1,4]
     BIC[l,4] <- BIC(fit4)
  }
}

######## p-values: individual model
p1_1 <- ifelse(p1<0.05,1,0)
apply(p1_1,2,mean)
p1_2 <- apply(p1_1,1,function(x){paste(x,collapse="_")})
length(unique(p1_2))

######## p-values: partialing model
p_a_1 <- ifelse(p_a<0.05,1,0)
apply(p_a_1,2,mean)
p_a_2 <- apply(p_a_1,1,function(x){paste(x,collapse="_")})
length(unique(p_a_2))

######## p-values: model with only composite M
p_b_1 <- ifelse(p_b<0.05,1,0)
apply(p_b_1,2,mean)

######## p-values: 1-level model
p2_1 <- ifelse(p2<0.05,1,0)
apply(p2_1,2,mean)
p2_2 <- apply(p2_1,1,function(x){paste(x,collapse="_")})
length(unique(p2_2))

######## p-values: 2-level model
p3_1 <- ifelse(p3<0.05,1,0)
apply(p3_1,2,function(x){mean(x,na.rm=T)})
p3_2 <- apply(p3_1,1,function(x){paste(x,collapse="_")})
length(unique(p3_2))


######## p-values: 4-level model
p4_r <- p4
for(l in 1:sample_number){
  p4_r[l,] <- ifelse(any(is.na(p4[l,])),rep(NA,3),p4[l,])
}

p4_r <- ifelse(p4_r<0.05,1,0)
table(p4_r[,1])
table(p4_r[,2])
table(p4_r[,3])


######### Model selection by BIC

choice <- apply(BIC,1,function(x){match(min(x,na.rm=T),x)})
table(choice)



