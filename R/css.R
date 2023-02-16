# standardize Matrix to have 0 mean and I var
stand <- function(x){
  n<-nrow(x)
  p<-ncol(x)
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- t(x1) %*% (x1)/(n-1)
  eva <- eigen(sigma)$values
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1/sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}

# power of. matrix
matpower <- function(a,alpha){
  small <- .00000001
  if (length(c(a))==1) {return(a^(alpha))} else {
    p1<-nrow(a)
    eva<-eigen(a)$values
    eve<-eigen(a)$vectors
    eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
    index<-(1:p1)[abs(eva)>small]
    evai<-eva
    evai[index]<-(eva[index])^(alpha)
    ai<-eve%*%diag(evai,length(evai))%*%t(eve)
    return(ai)}
}

# Calculate etaopg
## d = LL
## h = bandwidth in kernel mean estimate (default = 2)
opg <- function(x0,y,d,h){
  n <- nrow(x0)
  p <- ncol(x0)
  x <- stand(x0)
  b <- numeric(0)
  for(i in 1:n){
    del <- cbind(1,t(t(x) - x[i,]))
    w <- kern(x,h,i)
    bi <- wls(del,y,w)[2:(p+1)]
    b <- cbind(b,bi)
  }
  v1 <- eigen(b%*%t(b))$vectors[,1:d]
  v3 <- matpower(cov(x0),-0.5)%*%v1
  return(v3)
}

# compute kernel
kern <- function(x,h,i){
  del <- t(t(x)-x[i,])
  ndel <- diag(del%*%t(del))
  w <- (1/sqrt(2*pi))*(1/h)*exp((-1/2)*ndel/h^2)
  w <- c(w)
  w1 <- w/sum(w)
  return(w1)
}

# weighted least squares coefficients
wls <- function(x,y,w){
  b <- matpower(t(x*w)%*%(x),-1)%*%(t(x*w)%*%y)
  return(b)
}

# Calculate CSSDR
cssdr <- function(p,n,d,etaopg,x0,y,nslices){
  ## n=TT
  ## p=KK
  ## d=LL
  ## x0<-hFF
  ## y=yy
  ## objective function
  obj_cssdr <- function(veta){
    eta <- rbind(diag(1,d),matrix(veta,p-d,d))
    xc <- t(t(x0)-apply(x0,2,mean))
    yc <- y-mean(y)
    n <- nrow(xc)
    nwithin <- n/nslices
    g<- gs(xc,eta)
    tmp1<- t(xc)%*% g %*% matpower(t(g)%*%g,-1)
    egamma<- apply(tmp1 %*% t(g),1,mean)
    egaga<- tmp1 %*% t(g)%*% xc/n
    exx <- t(xc)%*%xc/n
    term1 <- 0
    term2 <- matrix(0,p,p)
    term3 <- matrix(0,p,p)
    terma1 <- matrix(0,p,p)
    terma2 <- 0
    terma3 <- matrix(0,p,p)
    for (i in 1:nslices){
      tmp <- order(yc)[(nwithin*(i-1)+1):(nwithin*i)]
      onex<-xc[tmp,]
      tmpg<- gs(onex,eta)
      tmpegamma <- apply(tmp1 %*% t(tmpg),1,mean)
      tmpegaga <- tmp1%*%t(tmpg)%*%tmpg%*%t(tmp1)/nwithin
      term1 <- term1+t(tmpegamma)%*% tmpegamma
      term2 <- term2+tmpegamma%*% t(tmpegamma)
      term3 <- term3+tmpegaga%*%tmpegaga
      tmpxx <- t(onex)%*%onex/nwithin
      tmpx <- apply(onex,2,mean)
      terma1 <- terma1+ tmpx%*% t(tmpegamma)
      terma2 <- terma2+ t(tmpegamma)%*% tmpx
      terma3 <- terma3+ tmpxx%*%tmpegaga
    }
    bb <- 2*(term1/nslices)*(term1/nslices)+sum(diag(2*(term2/nslices)%*%(term2/nslices)+2*term3/nslices-2*egaga%*%egaga))
    ab <- sum(diag(2*(terma1/nslices)%*%t(terma1/nslices)))+2*(terma2/nslices)*(terma2/nslices)+sum(diag(2*terma3/nslices-2*exx%*%egaga))
    obj<- bb-2*ab
    return (obj)
  }

  yc <- y-mean(y)
  xc <- t(t(x0)-apply(x0,2,mean))
  ### limit K ----
  if (d==1) {retg<-etaopg[-1]/etaopg[1]} else {
    retg<-etaopg[-c(1:d),]%*%matpower(etaopg[1:d,],-1)}
  outoptim1 <- optim(c(retg),obj_cssdr)
  reta <- matrix(outoptim1$par,p-d,d)
  eta <- rbind(diag(1,d),matrix(reta,p-d,d))
  return(eta)
}

# compute the g function in objective
gs <- function(xc,eta){
  d=ncol(eta)
  if (d==1) {g<-cbind(1,xc%*%c(eta),(xc%*%c(eta))^2)}
  if (d==2) {
    g <- cbind(1,xc%*%c(eta[,1]),
               xc%*%c(eta[,2]),
               (xc%*%c(eta[,1]))^2,
               (xc%*%c(eta[,1]))*(xc%*%c(eta[,2])),
               (xc%*%c(eta[,2]))^2)
  }
  return(g)
}

# predict
customGAM <- function(yy,XX){
  T0 <- length(yy)
  XX <- as.matrix(XX)
  TT <- nrow(XX)
  pp <- ncol(XX)
  XX <- as.data.frame(XX)
  colnames(XX) <- paste0("x",1:pp)
  xtemp <- as.data.frame(XX[1:T0,])
  colnames(xtemp) <- paste0("x",1:pp)
  gam.fit <- gam::gam(as.formula(paste0("yy~",paste0("s(x",1:pp,")",collapse="+"))),data=xtemp)
  preds <- predict(gam.fit,newdata=XX[TT,,drop=FALSE])
  return(preds[[1]])
}











