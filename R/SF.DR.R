#' Directional regression for sufficient forecasting
#'
#' @param y Response, T by 1 matrix
#' @param X Predictors, p by T matrix
#' @param newX New predictors, a vector contains p entries (or \code{NULL})
#' @param K The number of common factors (default = obtained
#' by \code{\link{getK}})
#' @param L The number of predictive indices, L is required to be no greater than
#' K (default = 1)
#' @param etaopg hyperparameter in DR (default = obtained by \code{opg})
#' @param nslices hyperparameter in DR (default = 3)
#'
#' @return Out-of-sample forecast for \code{newX}; or in-sample forecast for the last
#' observed data point if \code{newX} is \code{NULL}
#' @import gam
#' @export
#' @references
#' Luo, W., Xue, L., Yao, J. and Yu, X. (2022), Inverse moment methods for sufficient
#' forecasting using high-dimensional predictors, \emph{Biometrika} 109(2), 473–487.
#' @examples
#' utils::data(dataExample,package = "sufficientForecasting")
#' SF.DR(dataExample$y,dataExample$X,dataExample$newX)
SF.DR <- function(y, X, newX = NULL, K = "default", L = 1, etaopg = "default",
                  nslices = 3){
  # default K
  if(K == "default"){
    K <- getK(y, X, 12)
  }
  # warning
  ## format
  if(!is.matrix(X) | !is.matrix(y)){
    stop("X and y must be matrices")
  }
  ## X PP by TT
  ## y TT by 1
  if(dim(y)[1] != dim(X)[2]){
    stop("X must be a P by T matrix and y must be a T by 1 matrix")
  }
  ## newX, p by 1 vector
  pp <- nrow(X)
  if(!is.null(newX)){
    if(!is.vector(newX) | length(newX) != pp){
      stop("new predictors must be a vector containing p entries")
    }
  }
  ## L <= K & int & >= 1
  if(L > K | L < 1 | L%%1 != 0){
    stop("invalid L: L must be an integar and must be not smaller than 1 and
         not greater than K ")
  }
  ## K
  if(!is.numeric(K)){
    stop("invalid K: try K = 'default'")
  }
  ## ---SF.DR doesn't work when k = 1---
  if(K < 2 | K%%1 != 0){
    stop("invalid K: K must be an integar and not smaller than 2")
  }
  ## nslices
  maxi <- max(L,2)
  if(nslices < maxi | nslices%%1 != 0){
    stop("invalid nslices: nslices must be an intergar and >= max{L,2} is required")
  }
  ## etaopg
  if(length(etaopg) == 1 && etaopg != "default"){
    stop("invalid etaopg: try etaopg = 'default'")
  }
  if(length(etaopg) > 1){
      if(!is.matrix(etaopg)){
        stop("etaopg must be a K by L matrix")
  }
      if(dim(etaopg)[1] != K | dim(etaopg)[2] != L){
        stop("etaopg must be a K by L matrix")
  }
  }

  # in-sample forecasting
  if(is.null(newX)){
    p <- nrow(X)
    tt <- ncol(X)
    cy <- y[-length(y)]
    ny <- length(cy)

    PCA = eigen( t(X) %*% X )
    hFF = as.matrix(PCA$vectors[,1:K] * sqrt(tt))   # tt*KK
    hBB = X %*% hFF / tt
    ## default etaopg
    if(length(etaopg) == 1){
      etaopg <- opg(hFF[1:ny,,drop=FALSE], cy, L, h=2)
    }
    Phi.CSSDR = cssdr(p = K, n = tt,d = L,etaopg = etaopg,
                      x0 = hFF[1:ny,,drop=FALSE], y = cy, nslices)
    Predictor.CSS <- hFF %*% Phi.CSSDR
    res <- customGAM(cy,Predictor.CSS)
    return(round(res,4))
  }
  # out-sample forecasting
  if(!is.null(newX)){
    cX <- cbind(X,newX)
    ny <- length(y)
    p <- nrow(cX)
    tt <- ncol(cX)

    PCA = eigen( t(cX) %*% cX )
    hFF = as.matrix(PCA$vectors[,1:K] * sqrt(tt))   # tt*KK
    hBB = cX %*% hFF / tt
    ## default etaopg
    if(length(etaopg) == 1){
      etaopg <- opg(hFF[1:ny,,drop=FALSE], y, L, h=2)
    }
    Phi.CSSDR = cssdr(p = K, n = tt,d = L,etaopg = etaopg,
                      x0 = hFF[1:ny,,drop=FALSE], y = y, nslices)
    Predictor.CSS <- hFF %*% Phi.CSSDR
    res <- customGAM(y,Predictor.CSS)
    return(round(res,4))
  }






}
