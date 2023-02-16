# simulate data
DGP <- function(pp,TT,DGPpara){
  KK <- DGPpara$KK
  sig_UU <- DGPpara$sig_UU
  ## 1) Factors (auto corr)
  FF = matrix(0,TT,KK)
  FF.tr = rep(0,KK)
  for(tr in 1:50){FF.tr = 0.3*FF.tr+rnorm(KK)}
  FF[1,]=FF.tr
  for (iT in 2:TT){
    FF[iT,] = 0.3*FF[iT-1,] + 0.95*rnorm(KK)
  }
  ## 2) Loadings
  BB = matrix(0,pp,KK)
  for (iK in 1:KK){
    BB[,iK] = runif(pp,-1,2)
  }
  ## 3) Errors in the factor model
  UU = matrix(0,pp,TT)
  UU.tr = rep(0,pp)
  for(tr in 1:50)
  {UU.tr=0.5*UU.tr+sig_UU * rnorm(pp)}  # initialize time series
  UU[,1] = UU.tr
  for (iT in 2:TT){
    UU[,iT] = 0.5*UU[,iT-1] + sig_UU * rnorm(pp)
  }
  ## 4) Predictors
  XX = matrix(0,pp,TT)
  for (iT in 1:TT){
    XX[,iT] = BB %*% FF[iT,] + UU[,iT]
  }
  res <- list()
  res$FF <- FF
  res$BB <- BB
  res$XX <- XX
  res$UU <- UU
  return(res)
}

# slice average
slav = function(x,y,yunit){
  n = nrow(x)
  p = ncol(x)
  nslice = length(yunit)
  xgy = matrix(0,nslice,p)
  for(i in 1:nslice){
    xgy[i,] = apply(as.matrix(x[y==yunit[i],]),2,mean)}
  return(xgy)
}

# slice proportions
slprob = function(y,yunit){
  n = length(y)
  nslice = length(yunit)
  out = rep(0,nslice)
  for(i in 1:nslice){
    out[i] = length(y[y==yunit[i]])/n}
  return(out)
}

# discretize
discretize = function(y,yunit){
  n = length(y)
  y = y + .00001*mean(y)*rnorm(n)
  nsli = length(yunit)
  yord = y[order(y)]
  n = length(y)
  nwith = floor(n/nsli)
  divpt = rep(0,nsli-1)
  for(i in 1:(nsli-1)){
    divpt[i] = yord[i*nwith+1]}
  y1 = rep(0,n)
  y1[y >= divpt[nsli-1]] = nsli
  y1[y < divpt[1]] = 1
  for(i in 2:(nsli-1)){
    y1[(y >= divpt[i-1])&(y < divpt[i])] = i}
  return(y1)
}

#  sliced inverse regression
##  y must be discretized according to yunit for continuous y;
##  for discrete y, use original y with yunit being the distinct values in
##  discrete y
sir.cov <- function(x, y, discretization=T, nslice){
  if (discretization == T){
    yunit = 1:nslice
    ydis = discretize(y,yunit)
  }
  if (discretization == F){
    yunit = unique(y)
    ydis = y
  }

  xgy = slav(x,ydis,yunit)
  prob = slprob(ydis,yunit)
  p = ncol(x)
  out = matrix(0,p,p)
  nslice = length(yunit)
  for(i in 1:nslice){
    out = out + prob[i]*xgy[i,]%*%t(xgy[i,])}
  return(out)
}
