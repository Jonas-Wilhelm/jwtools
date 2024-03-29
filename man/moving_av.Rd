% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{moving_av}
\alias{moving_av}
\title{Calculate moving average}
\usage{
moving_av(x, window = 5)
}
\arguments{
\item{x}{Vector for which to calculate moving average.}

\item{window}{Width of the window for the moving average. Number of
neighboring values that will be averaged at each position of x.}
}
\description{
This function calculates a moving average of a vector for smoothing data
(e.g. time series).
}
\examples{
moving_av(x = rnorm(100), window = 5)
}
