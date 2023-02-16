##---code to prepare `testy` dataset goes here
# parameters
param <- matrix(data = c(1,100,101,2,6,1),nrow = 1,ncol = 6)
colnames(param) <- c("id","p","T","L","K","sigU")
param <- as.data.frame(param)
# Model1 simulate
set.seed(321)
LL = as.numeric(param["L"])
KK = as.numeric(param["K"])
DGPpara =list()
DGPpara$KK = KK
DGPpara$sig_UU = as.numeric(param["sigU"])
pp = as.numeric(param["p"])
TT = as.numeric(param["T"])
DGPdata <- DGP(pp,TT,DGPpara)
FF <- DGPdata$FF
XX <- DGPdata$XX
BB <- DGPdata$BB
yy = 0.8*FF[,1]+0.5*FF[,2]+0.3*FF[,3] + rnorm(TT)
y <- as.matrix(yy[1:100])

##---code to prepare `testX` dataset goes here
# parameters
param <- matrix(data = c(1,100,101,2,6,1),nrow = 1,ncol = 6)
colnames(param) <- c("id","p","T","L","K","sigU")
param <- as.data.frame(param)
# Model1 simulate
set.seed(321)
LL = as.numeric(param["L"])
KK = as.numeric(param["K"])
DGPpara =list()
DGPpara$KK = KK
DGPpara$sig_UU = as.numeric(param["sigU"])
pp = as.numeric(param["p"])
TT = as.numeric(param["T"])
DGPdata <- DGP(pp,TT,DGPpara)
XX <- DGPdata$XX
X <- XX[,1:100]

##---code to prepare `newX` dataset goes here
# parameters
param <- matrix(data = c(1,100,101,2,6,1),nrow = 1,ncol = 6)
colnames(param) <- c("id","p","T","L","K","sigU")
param <- as.data.frame(param)
# Model1 simulate
set.seed(321)
LL = as.numeric(param["L"])
KK = as.numeric(param["K"])
DGPpara =list()
DGPpara$KK = KK
DGPpara$sig_UU = as.numeric(param["sigU"])
pp = as.numeric(param["p"])
TT = as.numeric(param["T"])
DGPdata <- DGP(pp,TT,DGPpara)
XX <- DGPdata$XX
newX <- as.numeric(XX[,101])

##dataExample
dataExample <- list(X=X,y=y,newX=newX)
usethis::use_data(dataExample,overwrite = TRUE)




