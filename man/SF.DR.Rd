% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SF.DR.R
\name{SF.DR}
\alias{SF.DR}
\title{Directional regression for sufficient forecasting}
\usage{
SF.DR(y, X, newX = NULL, K = "default", L = 1, etaopg = "default", nslices = 3)
}
\arguments{
\item{y}{Response, T by 1 matrix}

\item{X}{Predictors, p by T matrix}

\item{newX}{New predictors, a vector contains p entries (or \code{NULL})}

\item{K}{The number of common factors (default = obtained
by \code{\link{getK}})}

\item{L}{The number of predictive indices, L is required to be no greater than
K (default = 1)}

\item{etaopg}{hyperparameter in DR (default = obtained by \code{opg})}

\item{nslices}{hyperparameter in DR (default = 3)}
}
\value{
Out-of-sample forecast for \code{newX}; or in-sample forecast for the last
observed data point if \code{newX} is \code{NULL}
}
\description{
Directional regression for sufficient forecasting
}
\examples{
utils::data(dataExample,package = "sufficientForecasting")
SF.DR(dataExample$y,dataExample$X,dataExample$newX)
}
\references{
Luo, W., Xue, L., Yao, J. and Yu, X. (2022), Inverse moment methods for sufficient
forecasting using high-dimensional predictors, \emph{Biometrika} 109(2), 473–487.
}
