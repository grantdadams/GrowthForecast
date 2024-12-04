#http://www.nielsensweb.org/fims2022/nn.R
load("nn.RData")
field <- tempMatFine
x <- as.numeric(rownames(field))
x1 <- x[row(field)]
y <- as.numeric(colnames(field))
x2 <- y[col(field)]
z <- as.vector(field)

library(TMB)
compile("nn.cpp")
dyn.load(dynlib("nn"))
data <- list()
data$noNodes <- 15
data$X <- cbind(x1,x2,1)
data$Y <- z

param <- list()
param$W1 <- matrix(0.1, nrow=ncol(data$X), ncol=data$noNodes)
param$W2 <- matrix(0.1, nrow=data$noNodes, ncol=1)

obj <- MakeADFun(data, param, DLL="nn", silent=TRUE)
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr,
                          control=list(eval.max=200000, iter.max=100000, trace=0)))

par(mfrow=c(1,3))
image(x, y, field, main="obs")
pred<-matrix(obj$report()$pred, nrow=nrow(field))
pred[is.na(field)]<-NA
image(x, y, pred, main="nn-field")

image(x, y, pred-field, main="obs and both as contours")
contour(x, y, pred, add=TRUE, lwd=2, col="darkgreen")
contour(x, y, field, add=TRUE, lwd=2, col="darkblue")

