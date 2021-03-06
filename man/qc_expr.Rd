% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{qc_expr}
\alias{qc_expr}
\title{A Function that Plots Quality Control for an Expression Data Matrix}
\usage{
qc_expr(data, USE_LOG2_FOR_EXPR = TRUE, ...)
}
\arguments{
\item{data}{A matrix of exrpression values for probe/gene (lines) and confditions (colums).}

\item{USE_LOG2_FOR_EXPR}{A boolean set to TRUE if we want to plot the log of the expression.}

\item{...}{Parameters passed to plot function.}
}
\description{
This function plots a boxplot for each sample in an expression data matrix.
}
