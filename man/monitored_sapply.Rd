% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{monitored_sapply}
\alias{monitored_sapply}
\title{A Monitored Version of `sapply`}
\usage{
monitored_sapply(vec, func, ...)
}
\arguments{
\item{vec}{A vector usualy gives as `sapply` first parameter.}

\item{func}{A function usualy gives as `sapply` second parameter.}

\item{...}{Parameters passed to `monitored_apply` function.}
}
\value{
A vector of the application of the `func` function to each element of the `vec` vector.
}
\description{
This function offer a monitored version of `sapply`.
}
