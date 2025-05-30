% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LapProbit.R
\name{lapord.fit}
\alias{lapord.fit}
\title{Approximate MLE of an ordinal data probit model with crossed random effects based on the
1st-order Laplace's approximation.}
\usage{
lapord.fit(x, y, f1, f2, zeta, gamma, lsA = log(0.5), lsB = log(0.5))
}
\arguments{
\item{x}{Model matrix.}

\item{y}{Response vector (binary data).}

\item{f1}{Factor representing the row effects.}

\item{f2}{Factor representing the columns effects.}

\item{zeta}{Numerical vector of estimated thresholds of an ordinal data probit model
for independent data.}

\item{gamma}{Numerical vector of estimated coefficients of an ordinal data probit
model for independent data.}

\item{lsA}{Numerical, starting value for logarithm of row effects standard deviation.}

\item{lsB}{Numerical, starting value for logarithm of column effects standard deviation.}
}
\value{
A list  containing the following components:
\item{\code{beta}}{the vector of estimated coefficients.}
\item{\code{alpha}}{the vector of estimated thresholds.}
\item{\code{sigmaA, sigmaB}}{the estimated standard deviation of
row and column random effects.}
}
\description{
Approximate MLE of an ordinal data probit model with crossed random effects based on the
1st-order Laplace's approximation.
}
\examples{
\dontrun{library(lme4)
mod.polr <- MASS::polr(factor(y) ~ studage + lectage + service + dept, method = "probit",
                       data = InstEval)
x <- model.matrix(mod.polr)[,-1]
fit.lap <- lapord.fit(x, InstEval$y, f1 = InstEval$s, f2 = InstEval$d, zeta = mod.polr$zeta,
                     gamma = mod.polr$coefficients)}

}
