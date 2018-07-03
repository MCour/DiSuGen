library(scales)


plotclustclin <- function(plotclin,datmix,clin,resflex){
  
  classrep <- table(apply(resflex@posterior$scaled,1,which.max))
  if(length(classrep)==1){
    Kuni <- as.numeric(names(classrep))
  }else{
    Kuni <- 1:length(classrep)
  }

  K <- resflex@k
  PAT <- unique(resflex@group)
  iPAT <- sapply(PAT,function(x){min(which(resflex@group==x))})
  clustres <- cbind(as.character(resflex@group[iPAT]),resflex@cluster[iPAT])
  transp2 <- 0.2
  
  test <- unlist(datmix[datmix$visit==2,names(datmix)==plotclin])
  nolong <- sum(is.na(test))==length(test)
  nbcomp <- length(table(datmix[,names(datmix)==plotclin]))
  compqt <- (nbcomp>10)
  compb <- (nbcomp==2)
  compp <- ((nbcomp<=10)&(nbcomp>2))
  
  par(mar=c(3.7,3.7,0.4,0.2),mfcol=c(1,1),mgp=c(2.5,1,0))
  
  if(nolong){
    boxplot(as.numeric(datmix[,names(datmix)==plotclin])~resflex@cluster,col=(Kuni)+1,
            ylab=plotclin,xlab="cluster")
  }else{
    xlabi <- "time"
    time2 <- seq(round(range(datmix$time),1)[1],round(range(datmix$time),1)[2],10) 
    TIME2 <- cbind(1,time2,time2^2,time2^3,time2^4,time2^5,time2^6)
    time <- datmix$time
    TIME <- cbind(1,time,time^2,time^3,time^4,time^5,time^6)
    rgpl <- range(as.numeric(gsub(",",".",datmix[,names(datmix)==plotclin])),na.rm=T)
    plot(x=time2,y=rep(NA,length(time2)),ylim=rgpl,type="l",col=(1:K)+1,ylab=plotclin,lwd=2,lty=1,
         xlab=xlabi)
    transp <- 0.4
    for(pat in unique(datmix$PATIENT)){
      patid <- which(datmix$PATIENT==pat)
      lines(as.numeric(gsub(",",".",datmix[patid,names(datmix)==plotclin]))~datmix[patid,3], col=alpha(resflex@cluster[patid]+1, transp), type="b")
    }
    if(compqt){
        i <- which(clin==plotclin)
        mu <- mean(as.numeric(datmix[,names(datmix)==plotclin]),na.rm=T)
        sigma <- sd(as.numeric(datmix[,names(datmix)==plotclin]),na.rm=T)
        if(length(clin)==1){
          nPR <- length(resflex@components[[1]][[1]]@parameters$coef)-1
          vpl <- TIME2[,1:(1+nPR)]%*%parameters(resflex)[1:(nPR+1),]
          sigpl <- parameters(resflex)[nPR+2,]*sigma
        }else{
          nPR <- length(resflex@components[[1]][[i]]@parameters$coef)-1
          vpl <- TIME2[,1:(1+nPR)]%*%parameters(resflex)[[i]][1:(nPR+1),] 
          sigpl <- parameters(resflex)[[i]][nPR+2,]*sigma
        }
        vpl2 <- vpl*sigma+mu
      for(k in Kuni){
        lines(x=time2,y=vpl2[,k],col=k+1,lwd=2,lty=1)
          polygon(x=c(time2,rev(time2)),y=c(vpl2[,k]-sigpl[k],rev(vpl2[,k]+sigpl[k])),
                  col=alpha(k+1, transp2),border=NA)
      }
    }
  } 
  
  legend("topright",legend=c(paste("cluster",1:K),"CI 68%"),col=c((1:K)+1,NA),
         pch=c(rep(15,K),NA),fill=c(rep(NA,K),alpha(1, transp2)),border=NA,bty="n")
}


plotclustgenet <- function(plotgenet,resflex,reslogitW){
  
  classrep <- table(apply(resflex@posterior$scaled,1,which.max))
  if(length(classrep)==1){
    Kuni <- as.numeric(names(classrep))
  }else{
    Kuni <- 1:length(classrep)
  }
  
  K <- resflex@k
  rownames(reslogitW) <- dimnames(resflex@concomitant@x)[[2]]
  toplot2 <- unique(which(reslogitW!=0,arr.ind=T)[,1])
  
  if(length(toplot2)==0){
    print('no concomitant variable was selected')
  }else{
    par(mar=c(6,4,1,1))
    if(length(toplot2)==1){
      resnew <- reslogitW[toplot2,]
      toplot <- resnew[plotgenet]
      names(toplot) <- rownames(reslogitW)[toplot2]
      plot(y=resnew,x=rep(1,K),xaxt="n",xlab=NA,
           ylab="coefficients of association with the clusters",
           col=c(2:(K+1)),cex=1,pch=1:K,lwd=1,
           ylim=c(min(0,min(resnew)),
                  max(0,max(resnew))))
    }else{
      resnew <- reslogitW[toplot2,]
      sortk <- order(resnew[,plotgenet])
      toplot <- resnew[sortk,plotgenet]
      matplot(resnew[sortk,],xaxt="n",xlab=NA,
              ylab="coefficients of association with the clusters",
              col=c(2:(K+1)),cex=1,pch=1:K,lwd=1,
              ylim=c(min(0,min(resnew)),
                     max(0,max(resnew))))
    }
    axis(1,labels=names(toplot),las=2,at=1:length(toplot),cex.axis=0.6)
    legend("topleft",legend=paste("cluster",1:K),col=c(2:(K+1)),pch=1:K,bty="n")
  } 
}
