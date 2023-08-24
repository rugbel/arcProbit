#' Laplace-based functions for categorical data.
#'
#' Approximate MLE of a binary data probit model with crossed random effects based on the
#' 1st-order Laplace's approximation.
#' @param x Model matrix.
#' @param y Response vector (binary data).
#' @param f1 Factor representing the row effects.
#' @param f2 Factor representing the columns effects.
#' @param gamma Numerical vector of estimated coefficients of binary data probit model.
#' @return A list  containing the following components:
#' \item{\code{beta}}{the vector of estimated coefficients.}
#' \item{\code{alpha}}{the vector of estimated thresholds.}
#' \item{\code{sigmaA, sigmaB}}{the estimated standard deviation of
#' row and column random effects.}
#' @export
#' @importFrom TMB MakeADFun
#' @importFrom ucminf ucminf
#' @examples
#' library(glmm)
#' data(salamander)
#' y <- salamander$Mate
#' mod.glm <- glm(Mate ~ Cross - 1, family = binomial("probit"), data = salamander)
#' x <- model.matrix(mod.glm)
#' fit.lap <- lapbin.fit(x, y, f1 = salamander$Male, f2 = salamander$Female, gamma = mod.glm$coef)
#'
lapbin.fit <- function(x, y, f1, f2, gamma, lsA = log(0.5), lsB = log(0.5)){
  g1 <- as.numeric(f1) - 1
  g2 <- as.numeric(f2) - 1
  p <- ncol(x)
	data.list <- list(X = x, y = y, g1 = g1, g2 = g2)
  parameters <- list(beta = gamma, lsigmaA = lsA, lsigmaB = lsB,
                     a = numeric(length(unique(g1))),
                     b = numeric(length(unique(g2))))
  obj <- TMB::MakeADFun(data = data.list, parameters = parameters,
                        DLL = "binCross", silent  = TRUE,  random = c("a", "b"))
  mle <- ucminf::ucminf(obj$par, obj$fn, obj$gr)
  betaAB <- mle$par[1:p]
	sigmaA <- exp(mle$par[p+1])
	sigmaB <- exp(mle$par[p+2])
  return(list(beta = betaAB, sigmaA = sigmaA, sigmaB = sigmaB))
}


#' Approximate MLE of an ordinal data probit model with crossed random effects based on the
#' 1st-order Laplace's approximation.
#' @param x Model matrix.
#' @param y Response vector (binary data).
#' @param f1 Factor representing the row effects.
#' @param f2 Factor representing the columns effects.
#' @param gamma Numerical vector of estimated coefficients of an ordinal data probit
#' model for independent data.
#' @param zeta Numerical vector of estimated thresholds of an ordinal data probit model
#' for independent data.
#' @return A list  containing the following components:
#' \item{\code{beta}}{the vector of estimated coefficients.}
#' \item{\code{alpha}}{the vector of estimated thresholds.}
#' \item{\code{sigmaA, sigmaB}}{the estimated standard deviation of
#' row and column random effects.}
#' @export
#' @importFrom TMB MakeADFun
#' @importFrom ucminf ucminf
#' @import lme4
#' @importFrom MASS polr
#' @examples
#' \dontrun{library(lme4)
#' mod.polr <- MASS::polr(factor(y) ~ studage + lectage + service + dept, method = "probit",
#'                        data = InstEval)
#' x <- model.matrix(mod.polr)[,-1]
#' fit.lap <- lapord.fit(x, InstEval$y, f1 = InstEval$s, f2 = InstEval$d, zeta = mod.polr$zeta,
#'                      gamma = mod.polr$coefficients)}
#'
lapord.fit <- function(x, y, f1, f2, zeta, gamma,  lsA = log(0.5), lsB = log(0.5)){
  g1 <- as.numeric(f1) - 1
  g2 <- as.numeric(f2) - 1
  p <- ncol(x)
  d <- length(unique(y))
  data.list <- list(X = x, y = y, d = d, maxtau = 100,
                    g1 = g1, g2 = g2)
  parameters <- list(tau = zeta, beta = gamma, lsigmaA = lsA, lsigmaB = lsB,
                     a = numeric(length(unique(g1))),
                     b = numeric(length(unique(g2))))
  obj <- TMB::MakeADFun(data = data.list, parameters = parameters,
                        DLL = "ordCross", silent  = TRUE,  random = c("a", "b"))
  mle <- ucminf::ucminf(obj$par, obj$fn, obj$gr)
  betaAB <- mle$par[d:(d+p-1)]
  alphaAB <- mle$par[1:(d-1)]
  sigmaA <- exp(mle$par[d+p])
  sigmaB <- exp(mle$par[d+p+1])
  return(list(beta = betaAB,  alpha = alphaAB,
              sigmaA = sigmaA, sigmaB = sigmaB))
}
