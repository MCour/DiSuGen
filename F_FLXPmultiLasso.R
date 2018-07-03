# A driver developped to add a Lasso selection 
# to the multinomial logistic regression for concomitant variables in the flexmix package
# see pages 25, 30 and 31 of 
# Grun, B., & Leisch, F. (2008). FlexMix version 2: finite mixtures with concomitant variables and varying and constant parameters.

# To test the functions :
# load("parameters.RData")
# x <- G
# y <- Z
# w <- rep(1,nrow(x))



FLXPmultiLasso <- function(formula = ~ 1) {
    z <- new("FLXPmultinom", name = "myConcomitant", formula = formula)
    multinom.fit2 <- function(x, y, w, ...) {
      r <- ncol(x)
      p <- ncol(y)
      if (p < 2) stop("Multinom requires at least two components.")
      # multinomial logistic regression with Lasso penalization:
      cv.glmnet(y=y,x=x,family="multinomial",alpha=1,nfolds=15,intercept=F)
      # one may choose a different alpha value to get a Lasso penalty mixed with a ridge penalty 
    }
    z@fit <- function(x, y, w, ...) {
            resglmnet <- multinom.fit2(x, y, w, ...)
            r <- ncol(x)
            p <- ncol(y)
            BETA <- matrix(data=NA,ncol=p,nrow=r)
            for(k in 1:p){
              BETA[,k] <-  as.matrix(coef(resglmnet, s = "lambda.min")[[k]])[-1]
            }
            AVG <- exp(x%*%BETA)
            pds <- AVG/rowSums(AVG)
            rownames(pds) <- rownames(x)
            return(pds)
    }
    z@refit <- function(x, y, w, ...) {
          resglmnet <- multinom.fit2(x, y, w, ...)
          r <- ncol(x)
          p <- ncol(y)
          BETA <- matrix(data=NA,ncol=p,nrow=r)
          for(k in 1:p){
            BETA[,k] <-  as.matrix(coef(resglmnet, s = "lambda.min")[[k]])[-1]
          }
          BETAp <- as.matrix(BETA[,-1]-BETA[,1])
          rownames(BETAp) <- colnames(x)
          betres <- lm(1~1) # only to get an object for which coef() is defined
          # (in order to use the FLXfillConcomitant flexmix method
          # for FLXPmultinom objects)
          betres$coefficients <- t(BETAp) 
          return(betres)
    }
    z
  }
