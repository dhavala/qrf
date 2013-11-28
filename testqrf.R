rm(list=ls(all=TRUE))
setwd('F:/work/qrf')
require(randomForest)
require(Hmisc)
require(hash)
#require(compiler)
require(quantregForest)

set.seed(123)
simu=0
if(simu)
{
  n <- 30
  p <- 10

  # input n x p
  x<-matrix(runif(n*p),n,p);
  b <- matrix(runif(p),p,1);b[c(2,3)]<-0
  y <- as.vector(x%*%b)
} else{
  
  data(airquality)
  ## remove observations with mising values
  airquality <- airquality[ !apply(is.na(airquality), 1,any), ]
  ## number of remining samples
  n <- nrow(airquality)
  ## divide into training and test data
  indextrain <- sample(1:n,round(0.6*n),replace=FALSE)
  x     <- as.matrix(airquality[ indextrain,2:6])
  y     <- as.vector(airquality[ indextrain,1])
  p <- ncol(x)
  n <- nrow(x)
}

# quant reg parameters (options)
probs <- c(0.1,0.5,0.9);
mtry <- 3
ntree <- 500
nodesize <- 3
nPerm <- 1

qrf1 <- qrf.fit(x,y=y, probs=c(0.1,0.5,0.9),ntree=ntree,importance=FALSE,nPerm=1)
qrf2 <- randomForest(x,y=y)
qrf3 <- quantregForest(x,y=y)

pr1 <- predict.qrf(qrf1,xold=x,xnew=x)
pr2 <- predict(qrf2,x,predict.all=TRUE,nodes=TRUE)
pr3 <- predict(qrf3,x)

y1<-pr1$predicted[,2]
y2<-pr2$aggregate
y3<-pr3[,2]  
pdf('all_prediction.pdf')
plot(y,y1,xlab='y',ylab='yhat',pch=16)
abline(lm(y1~y),lwd=2,col='black')
points(y,y2,col='red',pch=3)
abline(lm(y2~y),lwd=2,col='red')
points(y,y3,col='blue',pch=4)
abline(lm(y3~y),lwd=2,col='blue')
legend('topleft',legend=c('rf','qrf','qrf-R'),col=c('red','blue','black'),lwd=c(2,2,2))
dev.off()


pr1 <- predict.qrf(qrf1,xold=x)
pr2 <- predict(qrf2)
pr3 <- predict(qrf3)

y1<-pr1[,2]
y2<-pr2
y3<-pr3[,2]  
pdf('oob_prediction.pdf')
plot(y,y1,xlab='y',ylab='yhat',pch=16)
abline(lm(y1~y),lwd=2,col='black')
points(y,y2,col='red',pch=3)
abline(lm(y2~y),lwd=2,col='red')
points(y,y3,col='blue',pch=4)
abline(lm(y3~y),lwd=2,col='blue')
legend('topleft',legend=c('rf','qrf','qrf-R'),col=c('red','blue','black'),lwd=c(2,2,2))
dev.off()

qrf1a <- qrf.fit(x,y=y, probs=c(0.1,0.5,0.9),ntree=ntree,importance=TRUE,nPerm=1)
qrf2a <- randomForest(x,y=y, importance=TRUE)
print(qrf1a$importance)
print(qrf2a$importance)
print(qrf1a$importanceSD)
print(qrf2a$importanceSD)


# compare with linear regression estimates

lm.obj <- lm(y~x)
qrfX <- qrf.fit(x,y=y, probs=c(0.1,0.5,0.9),ntree=ntree,importance=FALSE,nPerm=1)
yhat <- qrfX$predicted
estBeta <-rep(NA,p)
for (ii in seq(1,p))
{
  x1 <-x;x1[,ii]=x[,ii]+1
  qrfX1 <-  qrf.fit(x1,y=y, probs=c(0.1,0.5,0.9),ntree=ntree,importance=FALSE,nPerm=1)
  yhat1 <- qrfX1$predicted
  lm1.obj <- lm(yhat1~yhat)
  estBeta[ii] <- coefficients(lm1.obj)[1]
}

qrfX <- qrf.fit(x,y=y, probs=c(0.5),ntree=ntree,importance=FALSE,nPerm=1)



