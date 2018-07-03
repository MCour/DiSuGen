# Function and inputs:
# clustering on variables listed in "clin", which families are listed in "fams" ("gaussian", "binomial" or "poisson"),
# * if "bingenet"=T, with logistic weights depending on variables listed in "xnam" (included an intercept)
# * with "nPR" polynomial components at most,
# * with "Kmax" clusters if "KKmax"=T, with "Kmax" clusters at most if "KKmax"=F.
# "long": vector indicating for each variable of "clin" wether the variable is longitudinal.
# "Lassosel": choose TRUE if a Lasso selection should be performed on concomitant variables
# "Constr1": choose TRUE if the identifiability constraint should be respected in the logistic regression
# "EMtype": choose "hard" for a CEM, "weighted" for an EM (flexmix "classify" option)

# Outputs:
# * "resclust": flexmix() output  
# * "clustres": summary of clustering results: probability of belonging to each cluster for each patient



DiSuGen <- function(datmix,fams=NA,long=NA,clin,xnam,Kmax=3,nPR=1,bingenet=T,KKmax=T,Lassosel=T,nrep=50,itermax=250,Constr1=T,EMtype="hard"){
  ## To test the function:
  # fams=NA
  # long=NA
  # Kmax=2
  # bingenet=T
  # KKmax=T
  # Lassosel=T
  # nrep=50
  # itermax=250
  # Constr1=T
  # EMtype="hard"

  if(is.na(fams)){fams<-rep("gaussian",length(clin))} 
  if(is.na(long)){long<-rep(T,length(clin))} 
  
  G <- as.matrix(datmix[datmix$visit==1,xnam])

  fmla <- as.formula(paste("~ ",paste(paste0("I(time^",1:nPR,")"), collapse= "+"),"|PATIENT")) 
  # regression function on time, squared time,...
  # "|PATIENT" : several visits per patient, all belonging to the same cluster
  
  # control parameters of flexmix():
  contrut <- list(verbose = 5, iter.max = itermax, classify=EMtype, minprior=0)

  modlist <- list()
  for(i in 1:length(clin)){
    if(long[i]==F){regf <- "1"}else{regf <- "."}#long[i]==F: non-longitudinal variables
    modlist[[i]] <- FLXMRglm(as.formula(paste(clin[i],"~",regf)),family=fams[i])
    if(fams[i]=="binomial"){
      modlist[[i]] <- FLXMRglm(as.formula(paste("cbind(",clin[i],",1-",clin[i],")~",regf)),family=fams[i])
    }
  }
  
  if(bingenet==T){
    if(Lassosel==F){
      conco <- FLXPmultinom(formula = as.formula(paste("~ ",paste(xnam, collapse= "+"),"-1")))
    }else{ 
      conco <- FLXPmultiLasso(formula = as.formula(paste("~ ",paste(xnam, collapse= "+"),"-1")))
    }
  }else{conco <- NULL}
  
  if(KKmax==F){Ksearch <- 1:Kmax}else{Ksearch <- Kmax}
    
  FLEXres <- list()
  BICres <- NA*numeric(length(Ksearch))# BIC results for each possible number of clusters
  for(k in Ksearch){
    FLEXinit <- list()# flexmix() output for each initialization
    BICinit <- rep(NA,nrep)# corresponding BIC
    print(paste0("K=",k))
    for(repini in 1:nrep){
      print(paste0("Initialisation ",repini,"/",nrep))
      FLEXinit[[repini]] <- try(flexmix(fmla, 
                                   data = datmix, 
                                   k = k,
                                   model= modlist, 
                                   concomitant = conco,
                                   control = contrut)) 
      resflexmix <- FLEXinit[[repini]]
      if(class(FLEXinit[[repini]])!="try-error"){
                I <- length(unique(datmix$PATIENT))
        if(sum(long)!=length(clin)){
          penPR <- sum(unlist(as.matrix(sapply(FLEXinit[[repini]]@components,function(x){lapply(x,function(x2){1+sum(x2@parameters$coef!=0)})}))[long,]))
        }else{
          penPR <- sum(unlist(as.matrix(sapply(FLEXinit[[repini]]@components,function(x){lapply(x,function(x2){1+sum(x2@parameters$coef!=0)})}))))
        }
       # BIC computation:
        BICinit[repini]  <- (-2*as.numeric(summary(FLEXinit[[repini]])@logLik)+
                              sum(long==F)*log(I)*(2*k)+ # non-longitudinal coefficients
                              log(nrow(datmix))*penPR) # longitudinal coefficients
                              # penclust*log(I)) # concomitant variables - no need since already a Lasso selection (with cross-validation) 
      }else{ 
        BICinit[repini]  <- NA
      }
    }
    if(KKmax==F){
      FLEXres[[k]] <- FLEXinit[[which.min(BICinit)]]
      BICres[k] <- min(BICinit,na.rm=T) 
    }
  }
  if(KKmax==F){
    K <- which.min(BICres)
    resflexmix <- FLEXres[[K]]
    BICmin <- min(BICres,na.rm=T)
  }else{
    K <- Kmax
    resflexmix <- FLEXinit[[which.min(BICinit)]]
    BICmin <- min(BICinit,na.rm=T)
  }
  
  PAT <- unique(resflexmix@group) # patient list
  iPAT <- sapply(PAT,function(x){min(which(resflexmix@group==x))})
  clustres <- as.matrix(resflexmix@posterior$scaled[iPAT,]) # probability for each patient to belong to each class
  colnames(clustres) <- paste("Cluster",1:K) 
  rownames(clustres) <- as.character(resflexmix@group[iPAT])
  reslogit <- assoclog(G,clustres,Constr1,SDcalc=T)
  
  return(list(resflexmix=resflexmix,clustres=clustres,BICmin=BICmin,reslogit=reslogit))
  
}