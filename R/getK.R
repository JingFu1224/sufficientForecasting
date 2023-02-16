#' Estimate the number of common factors K
#'
#' @param y Response, T by 1 vector
#' @param X Predictors, p by T matrix
#' @param Kmax A prescribed upper bound that possibly increases with p and T
#' (default = 12)
#'
#' @return Estimate of K
#' @references
#' Bai, J., and Ng, S. (2002), Determining the number of factors in approximate
#' factor models, \emph{Econometrica} 70(1), 191-221.
#'
#' Li, H., Li, Q. and Shi, Y. (2017), Determining the number of factors when the
#' number of factors can increase with sample size, \emph{Journal of Econometrics}
#' 197(1), 76–86.
#' @export
#'
getK <- function(y, X, Kmax = 12) {
  # error
  ## K
  if(Kmax < 2 | Kmax%%1 != 0){
    stop("invalid Kmax: Kmax must be an integar and not smaller than 2")
  }
  pp = nrow(X)
  TT = ncol(X)
  PCA = eigen( t(X) %*% X )
  IC <- function(kk)
  {
    hFF = PCA$vectors[,1:kk] * sqrt(TT)
    hBB = X %*% hFF / TT
    loss = (norm(X - hBB %*% t(hFF), type = 'F')^2)/(pp*TT)
    temp = (pp + TT)/(pp*TT)
    penn = kk*temp*log(1/temp)
    return(log(loss)+penn)
  }
  Klist = c(2:Kmax)
  Kop = Klist[which.min(sapply(Klist, IC))]
  return(Kop)
}
