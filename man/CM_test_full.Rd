% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistictsUsedInComparison.R
\name{CM_test_full}
\alias{CM_test_full}
\title{Computing value of Cramér von Mises test (\eqn{CM_n}) for full samples}
\source{
R. B. D’Agostino, M. A. Stephens. "Goodness-of-Fit Techniques". New York: Marcel Dekker (1986)
}
\usage{
CM_test_full(X_vec, theta_hat, lambda_hat, n)
}
\arguments{
\item{X_vec}{numeric, the values of observations from uncensored sample}

\item{theta_hat}{a estimator of parameter theta (shape parameter)}

\item{lambda_hat}{a estimator of parameter lambda (scale parameter)}

\item{n}{integer, amount of observations in data}
}
\value{
numeric,a value of test statistic \eqn{CM_n}
}
\description{
Classic Cramér von Mises test (\eqn{CM_n}) applied to data without censoring.
}
\keyword{internal}
