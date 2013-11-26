require(randomForest)
require(Hmisc)
require(hash)

# this function will be internal, to calculate check loss
ncmp_evalCheckLoss <- function(err,probs=0.5)
{
  # given x, calculate its the check loss
  return( mean( err*(probs-(err<0)) ) )
}

evalCheckLoss <- cmpfun(ncmp_evalCheckLoss)

ncp_importance.qrf <- function(object,xold,xnew=NULL,nPerm=1,probs=probs)
{
  
  if(class(object)!='qrf') stop('supplied object is not a quantile randomForest object')
  if(class(object$rfobject)!='randomForest') stop('supplied object is not a randomForest object')
  if(is.null(object$rfobject$forest) ) stop('randomForest does not have forest stored. try running with keep.forest=TRUE')
  
  y <- object$rfobject$y
  q <- length(probs)
  p <- ncol(xold)
  n <- nrow(xold)
  ntree <- object$rfobject$ntree
  imp.mn <- matrix(0.0,p,q)
  imp.sd <- matrix(0.0,p,q)
  
  nodes.orig <- nodes.perm <- attr(predict(object$rfobject,xold,predict.all=TRUE,nodes=TRUE),'nodes')
  inbag <- object$rfobject$inbag
  for (pp in seq(1,p)){
    
    xnew <- kronecker(xold,rep(1,nPerm))
    perm.ind <- matrix(replicate(nPerm,sample(n,n,replace=F)),ncol=1)
    # create nPerm blocks of x where pp-th predictor is permuted
    xnew[,pp] <- x[perm.ind,pp]
    # get the nodes by trickling down each tree for each p-variarte covariate
    nodes.perm <- attr(predict(object$rfobject,xnew,predict.all=TRUE,nodes=TRUE),'nodes')
    
    # nodes.orig (n x ntree) contains the leaf-nodes of orginal x
    # nodes.perm (nPerm*n x ntree) contains the leaf-nodes of x with p-th predictor permuted
    
    # need to calculate the check loss of oob predictions
    chkls <- matrix(0.0,ntree,q)
    
    for (tt in seq(1,ntree))
    {
      #create a hash for this tree
      # hash is just a nice way of looking at tapply list (key-value pair)
      h <- hash()
      clust.heads <- tapply(y,nodes.orig[,tt],quantile,probs=probs)
      #clust.heads <- tapply(y,nodes.orig[,tt],quantile,probs=rep(0.5,q))
      values(h,keys=names(clust.heads)) <- clust.heads
      
      # get residuals based on orginal oob predictors at oob samples
      ind.oob <- which(inbag[,tt]==0)
      y.oob <-matrix(y[ind.oob],ncol=1)
      yhat.oob  <- matrix(t(values(h,keys=make.keys(nodes.orig[ind.oob,tt]))),ncol=q)
      resid.oob <- -sweep(yhat.oob,1,y.oob,'-')
      
      # get residuals based on permuted predicters at oob samples
      ind.oob <- kronecker(ind.oob,rep(nPerm,1))
      y.oob <-y[ind.oob]
      yhat.perm  <- matrix(t(values(h,keys=make.keys(nodes.perm[ind.oob,tt]))),ncol=q)
      resid.perm <- -sweep(yhat.perm,1,y.oob,'-')
      rm(h)
      
      
      chkls.oob <- matrix(0.0,1,q)
      chkls.perm <- matrix(0.0,nPerm,q)
      
      for(qq in seq(1,q))
      {
        chkls.oob[qq] <- evalCheckLoss(resid.oob[,qq],probs[qq])
        resid.tmp <- matrix(resid.perm[,qq],ncol=nPerm)
        for(np in seq(1,nPerm))
        {
          chkls.perm[np,qq] <- evalCheckLoss(resid.tmp[,np],probs[qq])
        }
      }
      chkls[tt,] <- colMeans(chkls.perm)-chkls.oob
      
    }
    imp.mn[pp,] <- colMeans(chkls)
    imp.sd[pp,] <- apply(chkls,2,sd)
  }
  
  out<-list(imp.mn=imp.mn,imp.sd=imp.sd)
}

ncp_predict.qrf <- function(object,xold=x,xnew=NULL,probs=c(0.1,0.5,0.9))
{
  if(!is.null(object$predicted) & is.null(xnew)){
    return(object$predicted)
  }
     
  nodes.orig <- attr(predict(object$rfobject,xold,predict.all=TRUE,nodes=TRUE),'nodes')
  if(is.null(xnew))
  {
    inbag <- object$rfobject$inbag
    # all in-bag nodes will be assigned zero
    nodes.new <- nodes.orig*(1-inbag)
    
  } else  {
    nodes.new <- attr(predict(object$rfobject,xnew,predict.all=TRUE,nodes=TRUE),'nodes')
    
  }
  #nodes.orig <- as.integer(nodes.orig)
  #nodes.new <- as.integer(nodes.new)
  
  m <- nrow(nodes.new)
  n <- nrow(nodes.orig)
  wt <- matrix(0,n,m)
  ntree <- object$rfobject$ntree
  y <- object$rfobject$y
  all.zero <- matrix(0,1,m)
  for (tt in seq(1,ntree))
  {
    
    nodes.tt.new <- nodes.new[,tt]
    nodes.tt.orig <- nodes.orig[,tt]
    uniq.nodes.new <- unique(nodes.tt.new)
    wt.map <- outer(nodes.tt.orig,uniq.nodes.new,'==')
    wt.tt.mean <- colSums(wt.map)
    all.tt.zero <- (wt.tt.mean==0)
    
    wt.tt.mean <- wt.tt.mean + all.tt.zero
    wt.map2 <- sweep(wt.map,2,wt.tt.mean,'/')
        
    names(all.tt.zero) <- uniq.nodes.new
    colnames(wt.map2) <- uniq.nodes.new
    
    all.zero <- all.zero + all.tt.zero[paste(nodes.tt.new)]
    wt.net <- wt.map2[,paste(nodes.tt.new)]
    wt <- wt+wt.net
    
  }
  wt <- sweep(wt,1,(ntree-all.zero),'/')
  #wt <- sweep(wt,1,ntree,'/')
  
  getPred <- function(wt.vec,probs=probs,y=y)
  {
    pred <- wtd.quantile(y,wt.vec,normwt=TRUE,probs=probs)
    return(pred)
    
  }
  yhat <- apply(wt,2,getPred,probs=probs,y=y)
  yhat <- t(yhat)
  colnames(yhat) <- paste(probs)
  # return an m x q matrix where each column is the q-th quantile prediction
  out <- list(predicted=yhat,wt=wt)
  return(out)

}

require(compiler)
enableJIT(3)
predict.qrf <- cmpfun(ncp_predict.qrf)
importance.qrf <- cmpfun(ncp_importance.qrf)


# tweaked quantregForest interface
ncp_qrf.fit <-  function(x,y, probs=c(0.1,0.5,0.9),mtry= ceiling(ncol(x)/3),nodesize= 10,ntree= 1000,importance=FALSE,nPerm=1){
  
  ## Some checks
  
  if(! class(y) %in% c("numeric","integer") )
    stop(" y must be numeric ")
  
  if(is.null(nrow(x)) || is.null(ncol(x)))
    stop(" x contains no data ")
  
  if(length(unique(y))<=4)
    stop(" The response variable y contains less than 5 unique values! Quantile Regression assumes a continuous response variable. ")
  
  
  if(length(unique(y))<10)
    warning(" The response variable y contains less than 10 unique values! Quantile Regression assumes a continuous response variable.")
  
  
  if(mtry < 1 || mtry > ncol(x)){
    warning(" The value of mtry is too low or high! Has been reset to default value.")
    mtry <- max( floor(ncol(x)/3) ,1)
  }
  
  if( nrow(x) != length(y) )
    stop(" predictor variables and response variable must contain the same number of samples ")
  
  if (any(is.na(x))) stop("NA not permitted in predictors")
  if (any(is.na(y))) stop("NA not permitted in response")
  
  
  
  ## Check for categorial predictors with too many categories (copied from randomForest package)
  if (is.data.frame(x)) {
    ncat <- sapply(x, function(x) if(is.factor(x) && !is.ordered(x))
      length(levels(x)) else 1)
  } else {
    ncat <- 1
  }
  
  maxcat <- max(ncat)
  if (maxcat > 32)
    stop("Can not handle categorical predictors with more than 32 categories.")
  
  # check if the quantiles are ok
  if(! class(probs) %in% c("numeric"))  stop(" probs must be numeric ")
  if (any(probs<0)) stop("quantiles probs can not be negative")
  if (any(probs>1)) stop("quantiles probs can not be more than one")
  
  rf <- randomForest(x=x,y=y,keep.forest=TRUE, mtry=mtry, nodesize=nodesize, ntree=ntree, keep.inbag=TRUE)
  
  qrf <- list(rfobject=rf)
  
  qrf$importance <- NULL
  class(qrf) <- 'qrf'
  
  # compute vraible importance at prdeictors used in the training data
  if(importance)
  {
    importance.object <- importance.qrf(qrf,xold=x,nPerm=nPerm,probs=probs);
    qrf$importance <- importance.object$imp.mn
    qrf$importanceSD <- importance.object$imp.sd
  }
  
  # default compute oob predictions
  
  qrf$predicted <- NULL
  predicted.object <- predict.qrf(qrf,xold=x,probs=probs)
  qrf$predicted <- predicted.object$predicted
  qrf$wt <- predicted.object$wt
  
  return(qrf)
}

qrf.fit <- cmpfun(ncp_qrf.fit)

