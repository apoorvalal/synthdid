% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synthdid.R
\name{sc_estimate}
\alias{sc_estimate}
\title{synthdid_estimate for synthetic control estimates.
Takes all the same parameters, but by default, passes options to use the synthetic control estimator
By default, this uses only 'infinitesimal' ridge regularization when estimating the weights.}
\usage{
sc_estimate(Y, N0, T0, eta.omega = 0.000001, ...)
}
\arguments{
\item{Y}{the observation matrix.}

\item{N0}{the number of control units. Rows 1-N0 of Y correspond to the control units.}

\item{T0}{the number of pre-treatment time steps. Columns 1-T0 of Y correspond to pre-treatment time steps.}

\item{eta.omega}{determines the level of ridge regularization, zeta.omega = eta.omega * noise.level, as in synthdid_estimate.}

\item{...}{additional options for synthdid_estimate}
}
\value{
an object like that returned by synthdid_estimate
}
\description{
synthdid_estimate for synthetic control estimates.
Takes all the same parameters, but by default, passes options to use the synthetic control estimator
By default, this uses only 'infinitesimal' ridge regularization when estimating the weights.
}
