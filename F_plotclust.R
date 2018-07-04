library(scales)


plotclustclin <- function(clinvar,datmix,clin,resflex){
  
  classrep <- table(apply(resflex@posterior$scaled,1,which.max))
  if(length(classrep)==1){
    Kuni <- as.numeric(names(classrep))
  }else{
    Kuni <- 1:length(classrep)
  }

  Kmax <- resflex@k
  PAT <- unique(resflex@group)
  iPAT <- sapply(PAT,function(x){min(which(resflex@group==x))})
  clustres <- as.matrix(resflex@posterior$scaled[iPAT,]) # probability for each patient to belong to each class
  transp2 <- 0.2
  
  test <- unlist(datmix[datmix$visit==2,names(datmix)==clinvar])
  nolong <- sum(is.na(test))==length(test)
  nbcomp <- length(table(datmix[,names(datmix)==clinvar]))
  compqt <- (nbcomp>10)
  compb <- (nbcomp==2)
  compp <- ((nbcomp<=10)&(nbcomp>2))
  
  par(mar=c(3.7,3.7,0.4,0.2),mfcol=c(1,Kmax+1),mgp=c(2.5,1,0))
  
  if(nolong){
    boxplot(as.numeric(datmix[,names(datmix)==clinvar])~resflex@cluster,col=(Kuni)+1,
            ylab=clinvar,xlab="cluster")
  }else{
    xlabi <- "time"
    time2 <- seq(round(range(datmix$time),1)[1],round(range(datmix$time),1)[2],10) 
    TIME2 <- cbind(1,time2,time2^2,time2^3,time2^4,time2^5,time2^6)
    time <- datmix$time
    TIME <- cbind(1,time,time^2,time^3,time^4,time^5,time^6)
    rgpl <- range(as.numeric(gsub(",",".",datmix[,names(datmix)==clinvar])),na.rm=T)
    plot(x=time2,y=rep(NA,length(time2)),ylim=rgpl,type="l",col=(1:Kmax)+1,ylab=clinvar,lwd=2,lty=1,
         xlab=xlabi)
    transp <- 0.4
    if(compqt){
        i <- which(clin==clinvar)
        mu <- mean(as.numeric(datmix[,names(datmix)==clinvar]),na.rm=T)
        sigma <- sd(as.numeric(datmix[,names(datmix)==clinvar]),na.rm=T)
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
        lines(x=time2,y=vpl2[,k],col=k+1,lwd=2,lty=k)
          polygon(x=c(time2,rev(time2)),y=c(vpl2[,k]-sigpl[k],rev(vpl2[,k]+sigpl[k])),
                  col=alpha(k+1, transp2),border=NA)
      }
    }
    for(k in Kuni){
      par(mar=c(3.7,0.2,0.4,0.2),mgp=c(2.5,1,0))
      xlabi <- paste("cluster",k)
      xaxti <- "s"
      transp <- 0.3
      plot(x=time2,y=vpl2[,k],ylim=rgpl,type="l",col=k+1,
           ylab=NA,lwd=2,lty=k,cex.lab=1,cex.axis=1,
           xlab=xlabi,yaxt="n",xaxt=xaxti)
      for(pat in unique(datmix$PATIENT)){
        patid2 <- which(unique(datmix$PATIENT)==pat)
        patid <- which(datmix$PATIENT==pat)
        if(resflex@cluster[patid][1]==k){
          lines(as.numeric(gsub(",",".",datmix[patid,clinvar]))~datmix[patid,3], 
                col=alpha(Kuni[resflex@cluster[patid]]+1, clustres[patid2,k]*transp), 
                type="b") 
        }
      }
    }
  }
  legend("topleft",legend=c(paste("cluster",1:Kmax),"CI 68%"),col=c((1:Kmax)+1,NA),lty=c(1:Kmax,NA),
         pch=c(rep(NA,Kmax),NA),fill=c(rep(NA,Kmax),alpha(1, transp2)),border=NA,bty="n",cex=1)
}





plotclustgenet <- function(plotgenet,RESclust){
  
  k <- plotgenet
  resflex <- RESclust$resflexmix
  reslogitW <- RESclust$reslogit$West2
  toplot2 <- which(reslogitW[,k]!=0)
  resnew <- reslogitW[toplot2,]

  classrep <- table(apply(resflex@posterior$scaled,1,which.max))
  if(length(classrep)==1){
    Kuni <- as.numeric(names(classrep))
  }else{
    Kuni <- 1:length(classrep)
  }
  
  Kmax <- resflex@k
  # toplot2 <- unique(which(reslogitW!=0,arr.ind=T)[,1])
  
  if(length(toplot2)==0){
    print('no concomitant variable was selected')
  }else{
    par(mar=c(6,4,1,1))
    if(length(toplot2)==1){
      resnew <- reslogitW[toplot2,]
      toplot <- resnew[plotgenet]
      names(toplot) <- rownames(reslogitW)[toplot2]
      plot(y=resnew,x=rep(1,Kmax),xaxt="n",xlab=NA,
           ylab="coefficients of association with the clusters",
           col=c(2:(Kmax+1)),cex=1,pch=1:Kmax,lwd=1,
           ylim=c(min(0,min(resnew)),
                  max(0,max(resnew))))
    }else{
      par(bg="white",mgp=c(2.5,1,0),mar=c(3.4,9,0.1,0.1))
      resnew <- reslogitW[toplot2,]
      sortk <- order(resnew[,plotgenet])
      colpl <- (1:Kmax)+1
      toplot <- resnew[sortk,plotgenet]
      plot((1:length(sortk))~resnew[sortk,k],
           yaxt="n",ylab=NA,
           xlab=expression(omega * " coefficients"), 
           xaxt='s',
           col=colpl[k],cex=1,pch=k,lwd=1,
           xlim=c(min(0,min(resnew)),
                  max(0,max(resnew))))
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col ="snow1")
      abline(h=(1:length(sortk)),col="lightgray")
      abline(v=0,col="lightgray")
      for(kp in 1:Kmax){
        alphaplot <- 0.5
        if(k==plotgenet){alphaplot <- 1}
        points((1:length(sortk))~resnew[sortk,kp],col=alpha(colpl[kp],alphaplot),pch=kp)
      }
      axis(2,labels=names(toplot),cex.axis=0.8,
           las=2,at=1:length(toplot))
      
      legend("topleft",legend=paste("cluster",1:Kmax),col=colpl,pch=c(1:Kmax),
             bty="n",pt.cex=rep(1,Kmax))
      
      # SDmat <- RESclust$reslogit$SDmat[toplot2,]
      # sortk <- order(resnew[,k])
      # sd <- SDmat[sortk,k]*1.96
      # segments(y0=1:length(sortk),x0=resnew[sortk,k]-sd,x1=resnew[sortk,k]+sd,col=alpha(colpl[k],0.9))
      # epsilon <- 0.1
      # segments(y0=1:length(sortk)-epsilon,x0=resnew[sortk,k]-sd,y1=1:length(sortk)+epsilon,x1=resnew[sortk,k]-sd,col=colpl[k])
      # segments(y0=1:length(sortk)-epsilon,x0=resnew[sortk,k]+sd,y1=1:length(sortk)+epsilon,x1=resnew[sortk,k]+sd,col=colpl[k])
      # legend("topleft",legend=c(paste("cluster",1:Kmax),"95% CI"),col=c(colpl,1),pch=c(1:Kmax,NA),
      #        bty="n",pt.cex=c(rep(1,Kmax),NA),lty=c(rep(NA,Kmax),1))
    }

  } 
  
  
  
}