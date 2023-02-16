#' Conformal inference of the sufficient forecasting
#'
#' @param y Response, T by 1 matrix
#' @param X Predictors, p by T matrix
#' @param newX New predictors, a vector contains p entries (or \code{NULL})
#' @param type \code{LM} or \code{LLM} (default = \code{LM})
#' @param K The number of common factors (default = obtained
#' by \code{\link{getK}})
#' @param L The number of predictive indices, L is required to be no greater than
#' K (default = 1)
#' @param alpha Mis-coverage rate
#' @param discretization Hyperparameter in SIR (default = \code{TRUE})
#' @param nslices Hyperparameter in SIR (default = 10)
#'
#' @return A list with components
#' \describe{
#'   \item{yhat}{Out-of-sample forecast for \code{newX}; or in-sample forecast
#'   for the lastobserved data point if \code{newX} is \code{NULL}}
#'   \item{ci_lower}{Lower bound of conformal interval}
#'   \item{ci_upper}{Upper bound of conformal interval}
#' }
#' @export
#' @references
#' Yu, X., Yao, J. and Xue, L. (2022), Nonparametric estimation and conformal inference
#' of the sufficient forecasting with a diverging number of factors,
#' \emph{Journal of Business & Economic Statistics} 40(1), 342–354.
#' @examples
#' utils::data(dataExample,package = "sufficientForecasting")
#' SF.CI(dataExample$y,dataExample$X,type = "LM",alpha = 0.05)
SF.CI <- function(
    y, X, newX = NULL, type = "LM", K = "default", L = 1, alpha = 0.1,
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
  ## alpha
  if(alpha >= 1 | alpha <= 0){
    stop("invalid alpha value")
  }
  ## type error
  if(type != "LM" & type != "LLM"){
    stop("type must be one of 'LM' and 'LLM'")
  }
  ## alpha
  if(alpha >= 1 | alpha <= 0){
    stop("alpha must be (0,1)")
  }
  # pre
  SF_lm_fit = function(yy,Predictor){
    T0 = length(yy)
    Predictor = as.matrix(Predictor)
    xtemp = Predictor
    ytemp = yy
    beta = solve(t(xtemp)%*% xtemp)%*%(t(xtemp)%*%ytemp)
    eps = yy - Predictor %*% beta
    return(eps)
  }

  SF_LLR_fit = function(yy,Predictor){
    T0 = length(yy)
    Predictor = data.frame(Predictor)
    xtemp = Predictor
    ytemp = yy
    LLR.fit = loess(ytemp~.,data=xtemp,span=0.9, degree=1,
                    control=loess.control(surface = "direct"))
    eps = LLR.fit$residuals
    return(eps)
  }

  hateps <- function(yy, XX)
  {
    ## return the list of residuals
    ## input: yy: vector of length tt
    ##        XX: matrix of pp by tt

    ## PCA for factors and loadings
    tt = dim(XX)[2]
    PCA = eigen( t(XX) %*% XX )
    hFF = as.matrix(PCA$vectors[,1:K] * sqrt(tt))   # tt*KK
    hBB = XX %*% hFF / tt

    ## Predictors
    hFF.cov = sir.cov(hFF, yy, discretization, nslices)
    Phi.h = eigen(hFF.cov)$vectors[,1:L] # KK*LL
    Predictor = hFF %*% Phi.h

    ## return residuals
    eps_SF_LML = SF_lm_fit(yy,Predictor)
    eps_SF_LLRL = SF_LLR_fit(yy,Predictor)

    epsmat = cbind(eps_SF_LML, eps_SF_LLRL)
    colnames(epsmat) <- c("SF_LML", 'SF_LLRL')
    return(epsmat)
  }

  pyhat <- function(epsvec)
  {
    Tlast = length(epsvec)
    pvalue = mean(abs(epsvec) >= abs(epsvec)[Tlast])
    return(pvalue)
  }

  SF_conformal <- function(yy, XX)
  {
    Tlast = length(yy)
    y.grid = c(seq(min(yy)-2*sd(yy),max(yy)+2*sd(yy),length=200))

    # y-hat
    eps.hat = hateps(yy, XX)
    yhat = SF.SIR(y = y,X = X,newX = newX, type = type,
                  K = K, L = L, discretization = discretization,
                  nslices = nslices)

    # CI
    p.vec = matrix(NA,length(y.grid),2)
    colnames(p.vec) <- c('SF_LML', 'SF_LLRL')
    for (i in 1:length(y.grid)){
      yy[Tlast] = y.grid[i]
      eps.hat.mat = hateps(yy, XX)
      p.vec[i,] = apply(eps.hat.mat, 2, pyhat)
    }
    ci_SF_LML  = y.grid[p.vec[,'SF_LML']>alpha]
    ci_SF_LLRL = y.grid[p.vec[,'SF_LLRL']>alpha]

    ci_lower = c(min(ci_SF_LML), min(ci_SF_LLRL))
    ci_upper = c(max(ci_SF_LML), max(ci_SF_LLRL))
    resultmat = rbind(yhat, ci_lower, ci_upper)
    return(round(resultmat,4))
  }
  # in-sample
  if(is.null(newX)){
    ## LM
    out <- SF_conformal(y,X)
    if(type == "LM"){
      return(out[,1])
    }
    ## LLM
    if(type == "LLM"){
      return(out[,2])
    }
  }
  # out-of-sample
  if(!is.null(newX)){
    ## XX pp by tt+1
    ## y tt+1 by 1
    cX <- cbind(X,newX)
    cy <- c(y,mean(y))
    ## LM
    out <- SF_conformal(cy,cX)
    if(type == "LM"){
      return(out[,1])
    }
    ## LLM
    if(type == "LLM"){
      return(out[,2])
    }
  }
}
