#' Sliced inverse regression for sufficient forecasting
#'
#' @param y Response, T by 1 matrix
#' @param X Predictors, p by T matrix
#' @param newX New predictors, a vector contains p entries (or \code{NULL})
#' @param type \code{LM} or \code{LLM} (default = \code{LM}), \code{type = LM} fits
#' a linear regression of the response on the estimated predictive indices;
#' \code{type = LLM} fits a local linear regression
#' @param K The number of common factors (default = obtained
#' by \code{\link{getK}})
#' @param L The number of predictive indices, L is required to be no greater than
#' K (default = 1)
#' @param discretization Hyperparameter in SIR (default = \code{TRUE})
#' @param nslices Hyperparameter in SIR (default = 10)
#'
#' @return Out-of-sample forecast for \code{newX}; or in-sample forecast for the last
#' observed data point if \code{newX} is \code{NULL}
#' @import stats
#' @export
#' @references
#' Fan, J., Xue, L. and Yao, J. (2017), Sufficient forecasting using factor models,
#' \emph{Journal of econometrics} 201(2), 292–306.
#'
#' Yu, X., Yao, J. and Xue, L. (2022), Nonparametric estimation and conformal inference
#' of the sufficient forecasting with a diverging number of factors,
#' \emph{Journal of Business & Economic Statistics} 40(1), 342–354.
#' @examples
#' utils::data(dataExample,package = "sufficientForecasting")
#' SF.SIR(dataExample$y,dataExample$X,type = "LLM")
#'

SF.SIR <- function(y, X, newX = NULL, type = "LM", K = "default", L = 1,
                   discretization = TRUE, nslices = 10){
  # default K
  if(K == "default"){
    K <- getK(y, X, 12)
  }

  pp <- nrow(X)
  TT <- ncol(X)
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
  if(K < 1 | K%%1 != 0){
    stop("invalid K: K must be an integar and not smaller than 1")
  }
  ## nslices
  maxi <- max(L,2)
  if(nslices < maxi | nslices%%1 != 0){
    stop("invalid nslices: nslices must be an intergar and >= max{L,2} is required")
  }
  ## type error
  if(type != "LM" & type != "LLM"){
    stop("type must be one of 'LM' and 'LLM'")
  }
  ## discre
  if(!is.logical(discretization)){
    stop("discretization must be one of 'TRUE' and 'FALSE'")
  }
  # in-sample forecast
  if(is.null(newX)){
    ##  PCA for factors and loadings
    PCA = eigen( t(X) %*% X )
    hFF = as.matrix(PCA$vectors[,1:K] * sqrt(TT))   # tt*KK
    hBB = X %*% hFF / TT

    ##  Condition on hFF
    hFF.cov = sir.cov(as.matrix(hFF[-TT,]),y[-TT],discretization,nslices)
    Phi.h = eigen(hFF.cov)$vectors[,1:L] # KK*LL

    ##  Prediction
    Predictor = hFF %*% Phi.h    # tt*LL

    ## LM
    if(type == "LM"){
      SF_lm_forecasting = function(yy,Predictor){
        T0 = length(yy)
        Predictor = as.matrix(Predictor)
        xtemp = Predictor[1:(T0-1),]
        ytemp = yy[1:(T0-1)]
        beta = solve(t(xtemp)%*% xtemp)%*%(t(xtemp)%*%ytemp)
        hy = Predictor[T0,] %*% beta
        return(hy)
      }
      return(round(as.numeric(SF_lm_forecasting(y,Predictor)),4))
    }
    ## LLM
    if(type == "LLM"){
      SF_LLR_forecasting = function(yy,Predictor){
        T0 = length(yy)
        Predictor = data.frame(Predictor)
        xtemp = Predictor[1:(T0-1),]
        ytemp = yy[1:(T0-1)]
        LLR.fit = loess(ytemp~.,data=data.frame(xtemp),span=0.9, degree=1,
                        control=loess.control(surface = "direct"))
        hy = predict(LLR.fit,newdata=Predictor[T0,])[[1]]
        return(hy)
      }
      return(round(SF_LLR_forecasting(y,Predictor),4))
    }
  }
  # out-of-sample
  if(!is.null(newX)){
    ##  PCA for factors and loadings
    cX = cbind(X,newX)
    PCA = eigen( t(cX) %*% cX )
    hFF = as.matrix(PCA$vectors[,1:K] * sqrt(TT+1))   # (tt+1)*KK
    hBB = cX %*% hFF / (TT+1)

    ##  Condition on hFF
    hFF.cov = sir.cov(as.matrix(hFF[-(TT+1),]),y,discretization,nslices)
    Phi.h = eigen(hFF.cov)$vectors[,1:L] # KK*LL

    ##  Prediction
    Predictor = hFF %*% Phi.h    # tt*LL

    ## LM
    if(type == "LM"){
      SF_lm_forecasting = function(yy,Predictor){
        T0 = length(yy)
        Predictor = as.matrix(Predictor)
        xtemp = Predictor[1:T0,]
        ytemp = yy[1:T0]
        beta = solve(t(xtemp)%*% xtemp)%*%(t(xtemp)%*%ytemp)
        hy = Predictor[(T0+1),] %*% beta
        return(hy)
      }
      return(round(as.numeric(SF_lm_forecasting(y,Predictor)),4))
    }
    ## LLM
    if(type == "LLM"){
      SF_LLR_forecasting = function(yy,Predictor){
        T0 = length(yy)
        Predictor = data.frame(Predictor)
        xtemp = Predictor[1:T0,]
        ytemp = yy[1:T0]
        LLR.fit = loess(ytemp~.,data=data.frame(xtemp),span=0.9, degree=1,
                        control=loess.control(surface = "direct"))
        hy = predict(LLR.fit,newdata=Predictor[(T0+1),])[[1]]
        return(hy)
      }
      return(round(SF_LLR_forecasting(y,Predictor),4))
    }
  }
}
