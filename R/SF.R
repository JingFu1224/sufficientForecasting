#' Select a method from PC, SIR and DR to do point prediction
#'
#' @param y Response, T by 1 matrix
#' @param X Predictors, p by T matrix
#' @param newX New predictors, a vector contains p entries (or \code{NULL})
#' @param K The number of common factors (default = obtained by
#' \code{\link{getK}})
#' @param L The number of predictive indices, L is required to be no greater than
#' K (default = 1)
#' @param method Select one from \code{PC}, \code{SIR.LM}, \code{SIR.LLM}
#' and \code{DR} to do point prediction (default = \code{SIR.LM})
#' @param hyperparameter A list of parameters for the corresponding method
#' \describe{
#'   \item{For \code{SIR.LM} and \code{SIR.LLM}:}{arguments \code{discretization}
#'   and \code{nslices}. See \code{\link{SF.SIR}} for detail}
#'   \item{For \code{DR}:}{arguments \code{etaopg} and \code{nslices}.
#'   See \code{\link{SF.DR}} for detail}
#' }
#'
#' @return Out-of-sample forecast for \code{newX}; or in-sample forecast for the last
#' observed data point if \code{newX} is \code{NULL}
#' @references
#' Fan, J., Xue, L. and Yao, J. (2017), Sufficient forecasting using factor models,
#' \emph{Journal of econometrics} 201(2), 292–306
#'
#' Luo, W., Xue, L., Yao, J. and Yu, X. (2022), Inverse moment methods for sufficient
#' forecasting using high-dimensional predictors, \emph{Biometrika} 109(2), 473–487.
#'
#' Yu, X., Yao, J. and Xue, L. (2022), Nonparametric estimation and conformal inference
#' of the sufficient forecasting with a diverging number of factors,
#' \emph{Journal of Business & Economic Statistics} 40(1), 342–354.
#' @export
#' @examples
#' utils::data(dataExample,package = "sufficientForecasting")
#' SF(dataExample$y,dataExample$X,method = "SIR.LLM",
#' hyperparameter = list(nslices = 5,discretization = TRUE))
#' \dontrun{
#' SF(dataExample$y,dataExample$X,dataExample$newX,method = "DR")
#' SF(dataExample$y,dataExample$X,dataExample$newX,method = "PC")
#' }
SF <- function(y, X, newX = NULL, K = "default", L = 1,
                        method = "SIR.LM", hyperparameter = list()){
  # hyperparameters
  discretization <- hyperparameter$discretization
  nslices <- hyperparameter$nslices
  etaopg <- hyperparameter$etaopg

  # default discretization
  if(is.null(discretization)){
    discretization <- TRUE
  }
  # default nslices
  if(is.null(nslices)){
    if(method == "PC" || method == "SIR.LM" || method == "SIR.LLM"){
      nslices <- 10
    }
    if(method == "DR"){
      nslices <- 3
    }
  }
  # etaopg
  if(is.null(etaopg)){
    etaopg <- "default"
  }
  # error
  ## method
  if(method != "PC" & method != "SIR.LM" & method != "DR" & method != "SIR.LLM"){
    stop("'method' must be one of 'PC', 'SIR.LM', 'SIR.LLM' and 'DR' ")
  }

  # PC
  if (method == "PC"){
    res <- SF.PC(y,X,newX,K,L)
  }
  # SIR.LM
  if (method == "SIR.LM"){
    res <- SF.SIR(y,X,newX,type = "LM",K,L,discretization,nslices)
  }
  # SIR.LLM
  if (method == "SIR.LLM"){
    res <- SF.SIR(y,X,newX,type = "LLM",K,L,discretization,nslices)
  }
  # DR
  if (method == "DR"){
    res <- SF.DR(y,X,newX,K,L,etaopg,nslices)
  }
  return(res)
}
