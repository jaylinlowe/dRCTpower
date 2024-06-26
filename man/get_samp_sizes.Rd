% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_samp_sizes.R
\name{get_samp_sizes}
\alias{get_samp_sizes}
\title{Calculate sample sizes by subgroup}
\usage{
get_samp_sizes(Y, X, grouping_col, preds, effect_size, alpha, beta)
}
\arguments{
\item{Y}{a vector of responses}

\item{X}{matrix of covariates}

\item{grouping_col}{name of a column in X that defines subgroups}

\item{preds}{vector of predictions}

\item{effect_size}{desired minimal detectable difference}

\item{alpha}{desired Type II error rate}

\item{beta}{desired Type II error rate}
}
\value{
A data frame with 7 columns: \code{def}, \code{num}, \code{MSE}, \code{variance},\code{resid_var},
\code{samp_size}, and \code{samp_size_without}.
\code{def} is the definition of the subgroup, which may be an interval or a single
factor level, depending on the values in \code{grouping_col}. \code{num} is the number of
observations in that subgroup. \code{MSE} is the mean squared error of the predictions
for that subgroup. \code{variance} is the variance Y of the subgroup. \code{resid_var} is
the variance of the residuals calculated from \code{preds} and \code{Y}. \code{resid_var} is assumed
to be \eqn{\sigma^2} and is passed through the \code{samp_size} function. \code{samp_size}
is the resulting estimated sample size needed if the RCT population resembles this
subgroup. \code{samp_size_without} is the estimated sample size needed if the RCT
population resembles this subgroup and auxiliary data is not used to improve precision.
}
\description{
Calculate sample sizes by subgroup
}
