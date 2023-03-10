% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataExample.R
\docType{data}
\name{dataExample}
\alias{dataExample}
\title{A simulated dataset}
\format{
\itemize{
\item \code{dataExample$X} simulated predictors; a p by T matrix
\item \code{dataExample$y} simulated responses; a T by 1 vector
\item \code{dataExample$newX} simulated new predictors for forecasting; a p by 1 vector
}
}
\usage{
dataExample
}
\description{
This is a simulated dataset consisting of predictors, responses, and additional
predictors to make forecasting at a future time. In this example, the dimension
of predictors is p=100, and the number of observed time periods is T=100.
}
\keyword{datasets}
