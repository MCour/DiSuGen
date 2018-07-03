assoclog <- function(G,clustres,Constr1=T,SDcalc=T){
  
  K <- ncol(clustres)
  J <- ncol(G)
  BETA <- matrix(data=0,ncol=K,nrow=J)
  rownames(BETA) <- colnames(G)
  colnames(BETA) <- paste("Cluster",1:K)

  resglmnet <- cv.glmnet(y=clustres,x=G,family="multinomial",alpha=1,nfolds=15,intercept=F)
  lam <- which(resglmnet$lambda==resglmnet$lambda.min)
  for(k in 1:K){
    BETA[,k] <- as.matrix(resglmnet$glmnet.fit$beta[[k]])[,lam]
  }
  if(K==2|Constr1==T){BETA <- as.matrix(BETA-BETA[,1])}
  elim <- which(rowSums(BETA!=0)==K)
  for(i in elim){
    toelim <- order(abs(BETA[i,]))[c(1,2)]
    BETA[i,toelim] <- 0
  }

  clustresSNPs <- as.matrix(BETA)
  snpnn <- which(rowSums(abs(clustresSNPs))!=0)
  nbsnp <- length(snpnn)
  PVALS <- NA
  SDmat <- NA
  G2 <- G[,snpnn]
  clustresSNPs2 <- clustresSNPs[snpnn,]
  if(nbsnp!=0){
    
    if(nbsnp==1){
      mask <- as.vector(rbind(0,t(clustresSNPs2!=0)*1))
    }else{
      mask <- as.vector(rbind(0,(clustresSNPs2!=0)*1))
    }
    resnnet <- nnet(x=G2, y=clustres, size=0, mask = mask, skip=TRUE,
                    softmax=TRUE, Hess=T, rang=0)
    thetaest <- coef(resnnet)
    clustresSNPs2 <- matrix(data=thetaest,ncol=K)[-1,]
    if(nbsnp==1){
      clustresSNPs2 <- as.matrix(t(clustresSNPs2))
    }
    clustresSNPs[snpnn,] <- clustresSNPs2
    
    if(SDcalc==T){
      Hess <- resnnet$Hessian
      Hess2 <- Hess[mask==1,mask==1]
      InvHess <- try(solve(Hess2))
      if(class(InvHess)=="try-error"){
        PVALS <- NA
        SDmat <- NA
      }else{
        IncEst <- sqrt(diag(InvHess))
        IncEst2 <- thetaest*NA
        IncEst2[thetaest!=0] <- IncEst
        IncEst2[thetaest==0] <- 1 
        pvals <- pnorm(-abs(thetaest)/IncEst2,mean=0,sd=1)
        pvals[thetaest==0] <- NA
        IncEst2[thetaest==0] <- NA
        PVALS <- clustresSNPs*NA
        SDmat <- clustresSNPs*NA
        PVALS[snpnn,]<- matrix(data=pvals,ncol=ncol(clustresSNPs2))[-1,]
        SDmat[snpnn,]<- matrix(data=IncEst2,ncol=ncol(clustresSNPs2))[-1,]
        rownames(PVALS) <- colnames(G)
        rownames(SDmat) <- colnames(G)
        if(nbsnp>1){
          colnames(SDmat) <-  colnames(PVALS) <- colnames(BETA)
        }
      }
    }
    
    
  }
  
  return(list(West2=clustresSNPs,PVALS=PVALS,SDmat=SDmat,snpnn=snpnn))
}