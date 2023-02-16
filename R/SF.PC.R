#' Principal component regression for sufficient forecasting
#'
#' @param y Response, T by 1 matrix
#' @param X Predictors, p by T matrix
#' @param newX New predictors, a vector contains p entries (or \code{NULL})
#' @param K The number of common factors (default = obtained
#' by \code{\link{getK}})
#' @param L The number of principal components used in the prediction,
#' L is required to be no greater than K (default = \code{K})
#'
#' @return Out-of-sample forecast for \code{newX}; or in-sample forecast for the last
#' observed data point if \code{newX} is \code{NULL}
#' @export
#'
#' @examples
#' utils::data(dataExample,package = "sufficientForecasting")
#' SF.PC(dataExample$y,dataExample$X)
SF.PC <- function(y, X, newX = NULL, K = "default", L = "default"){
  # Estimate K\L
  if(K == "default"){
    K <- getK(y, X, 12)
  }
  if(L == "default"){
    L <- K
  }
  pp <- nrow(X)
  TT <- ncol(X)
  # error
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
  if(!is.null(newX)){
    if(!is.vector(newX) | length(newX) != pp){
      stop("new predictors must be a vector containing p entries")
    }
  }
  ## L <= K & int & >= 1
  if(!is.numeric(L)){
    stop("invalid L: try L = 'default'")
  }
  if(L > K | L < 1 | L%%1 != 0){
    stop("invalid L: L must be an integar and must be not smaller than 1 and
         not greater than K ")
  }
  ## K
  if(!is.numeric(K)){
    stop("invalid K: try K = 'default'")
  }
  if(K < 1 | K%%1 != 0){
    stop("invalid K: K must be an integar and not smaller than 1")
  }
  # in-sample forecasting
  if(is.null(newX)){
    SF_lm_forecasting = function(yy,Predictor){
      # one-step forecasting
      T0 = length(yy)
      Predictor = as.matrix(Predictor)
      xtemp = Predictor[1:(T0-1),]
      ytemp = yy[1:(T0-1)]
      beta = solve(t(xtemp)%*% xtemp)%*%(t(xtemp)%*%ytemp)
      hy = Predictor[T0,] %*% beta
      return(hy)
    }
    ##  PCA
    PCA = eigen( t(X) %*% X )
    hFF = as.matrix(PCA$vectors[,1:K] * sqrt(TT))   # tt*KK
    hBB = X %*% hFF / TT

    ##  Prediction
    return(round(as.numeric(SF_lm_forecasting(y, as.matrix(hFF[,1:L]))),4))
  }
  # out-of-sample forecasting
  if(!is.null(newX)){
    SF_lm_forecasting = function(yy,Predictor){
      # one-step forecasting
      T0 = length(yy)
      Predictor = as.matrix(Predictor)
      xtemp = Predictor[1:T0,]
      ytemp = yy
      beta = solve(t(xtemp)%*% xtemp)%*%(t(xtemp)%*%ytemp)
      hy = Predictor[(T0+1),] %*% beta
      return(hy)
    }
    ## PCA
    cX = cbind(X,newX)
    PCA = eigen( t(cX) %*% cX )
    hFF = as.matrix(PCA$vectors[,1:K] * sqrt(TT+1))   # tt+1*KK
    hBB = cX %*% hFF / (TT+1)

    ## Prediction
    return(round(as.numeric(SF_lm_forecasting(y, as.matrix(hFF[,1:L]))),4))
  }
}
