#' Functions for Pairwise Likelihood estimation for categorical probit data model.
#'
#' Fit a binary probit model with crossed random effects by the PL method.
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
#' @import Rcpp
#' @import RcppArmadillo
#' @import mvtnorm
#' @examples
#' library(glmm)
#' data(salamander)
#' y <- salamander$Mate
#' mod.glm <- glm(Mate ~ Cross - 1, family = binomial("probit"), data = salamander)
#' x <- model.matrix(mod.glm)
#' fit.pl <- plbin.fit(x, y, f1 = salamander$Male, f2 = salamander$Female,
#'                       gamma = mod.glm$coef)
#'
plbin.fit <- function(x, y, f1, f2, gamma){
  eta <- as.vector(x %*% gamma)
  ## now the two lists
  list1 <- split(1:nrow(x), f1)
  len1 <- unlist(lapply(list1, length))
  which.list1 <- which(len1 > 1)
  list1 <- list1[which.list1]
  list2 <- split(1:nrow(x), f2)
  len2 <- unlist(lapply(list2, length))
  which.list2 <- which(len2 > 1)
  list2 <- list2[which.list2]

  ## obtain the tetrachoric correlations
  rhoa <- optimise(lik2_binary,  lower = 0, upper = 1,  etaG = eta, list_ind = list1,
                   y = y, maximum = TRUE)$maximum
  rhob <- optimise(lik2_binary,  lower = 0, upper = 1,  etaG = eta, list_ind = list2,
                   y = y, maximum = TRUE)$maximum
  ## back to the conditional parameterization
  sigma2A <- rhoa / (1 - rhoa - rhob)
  sigma2B <- rhob / (1 - rhoa - rhob)
  betaAB <-  gamma * sqrt(1 + sigma2A + sigma2B)
  return(list(beta = betaAB, sigmaA = sqrt(sigma2A), sigmaB = sqrt(sigma2B)))
}


#' Fit an ordinal probit model with crossed random effects by the PL method.
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
#' @import Rcpp
#' @import lme4
#' @importFrom MASS polr
#' @import RcppArmadillo
#' @import mvtnorm
#' @examples
#' \dontrun{library(lme4)
#' mod.polr <- MASS::polr(factor(y) ~ studage + lectage + service + dept, method = "probit",
#'                        data = InstEval)
#' x <- model.matrix(mod.polr)[,-1]
#' fit.pl <- plord.fit(x, InstEval$y, f1 = InstEval$s, f2 = InstEval$d, zeta = mod.polr$zeta,
#'                    gamma = mod.polr$coefficients)}
#'
plord.fit <- function(x, y, f1, f2, zeta, gamma){
  eta <- as.vector(x %*% gamma)
  alpha <- c(-100, zeta, 100)
  ## now the two lists
  list1 <- split(1:nrow(x), f1)
  len1 <- unlist(lapply(list1, length))
  which.list1 <- which(len1 > 1)
  list1 <- list1[which.list1]
  list2 <- split(1:nrow(x), f2)
  len2 <- unlist(lapply(list2, length))
  which.list2 <- which(len2 > 1)
  list2 <- list2[which.list2]
  ## obtain the tetrachoric correlations
  rhoa <- optimise(lik2_ord,  lower = 0, upper = 1,  tauG = alpha, etaG = eta,
                   list_ind = list1, y = y, maximum = TRUE)$maximum
  rhob <- optimise(lik2_ord,  lower = 0, upper = 1,  tauG = alpha, etaG = eta,
                   list_ind = list2, y = y, maximum = TRUE)$maximum
   ## back to the conditional parameterization
  sigma2A <- rhoa / (1 - rhoa - rhob)
  sigma2B <- rhob / (1 - rhoa - rhob)
  betaAB <-  gamma * sqrt(1 + sigma2A + sigma2B)
  alphaAB <- alpha * sqrt(1 + sigma2A + sigma2B)
  return(list(beta = betaAB, alpha = alphaAB,
              sigmaA = sqrt(sigma2A),
              sigmaB = sqrt(sigma2B)))
}
