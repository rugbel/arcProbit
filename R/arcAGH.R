#' Fitting functions for categorical data.
#'
#' Fit a binary probit model with crossed random effects by the ARC method.
#' @param x Model matrix.
#' @param y Response vector (binary data).
#' @param f1 Factor representing the row effects.
#' @param f2 Factor representing the columns effects.
#' @param obj_glm Object representing the glm fit of a probit model
#' as returned by \code{glm}.
#' @param nq Number of quadrature points. Default is 25.
#' @param niter Number of Newton-Raphson iteration to locate
#' the mode of each random effect. Default is 10.
#' @param get_effects Returns the random effects prediction.
#' Default is FALSE.
#' @param get_se Returns the standard errors of beta coefficients.
#' Default is TRUE.
#' @return A list  containing the following components:
#' \item{\code{beta}}{the vector of estimated coefficients.}
#' \item{\code{beta.se}}{the estimated standard errors of the
#' coefficients.}
#' \item{\code{sigmaA, sigmaB}}{the estimated standard deviation of
#' row and column random effects.}
#' \item{\code{a.est, b.est}}{the predicted random effects.}
#' @export
#' @import glmm
#' @import Rcpp
#' @importFrom sandwich vcovCL
#' @importFrom statmod gauss.quad
#' @examples
#' library(glmm)
#' data(salamander)
#' y <- salamander$Mate
#' mod.glm <- glm(Mate ~ Cross - 1, family = binomial("probit"), data = salamander)
#' x <- model.matrix(mod.glm)
#' fit.arc <- arcbin.fit(x, y, f1 = salamander$Male, f2 = salamander$Female,
#'                       obj_glm = mod.glm)
#'
arcbin.fit <- function(x, y, f1, f2, obj_glm, nq = 25, niter = 10,
                       get_effects = FALSE, get_se = TRUE){
  ## transform ys to {-1, 1}
  w <- 2 * y - 1
  gamma <- obj_glm$coefficients
  eta <- as.vector(x %*% gamma)

  ## get the lists for taus
  listw1 <- split(w, f1)
  len1 <- unlist(lapply(listw1, length))
  sel1 <- which(len1 == 1)
  which.list1 <- which(len1 > 1)
  listw1 <- listw1[which.list1]
  listeta1 <- split(eta, f1)[which.list1]
  listw2 <- split(w, f2)
  len2 <- unlist(lapply(listw2, length))
  sel2 <- which(len2 == 1)
  which.list2 <- which(len2 > 1)
  listw2 <- listw2[which.list2]
  listeta2 <- split(eta, f2)[which.list2]

  ## obtain the tetrachoric correlations
  obj.gh <- statmod::gauss.quad(nq, "hermite")
  ws <- obj.gh$weights * exp(obj.gh$nodes^2)
  rhoa <- stats::optimise(likAGH, interval = c(0, 1), list_eta = listeta1,
                   list_w = listw1, niter = niter,
                   ws = ws, z = obj.gh$nodes, maximum = TRUE)$maximum
  rhob <- stats::optimise(likAGH, interval = c(0, 1), list_eta = listeta2,
                   list_w = listw2, niter = niter,
                   ws = ws, z = obj.gh$nodes, maximum = TRUE)$maximum

  ## back to the conditional parameterization
  sigma2A <- rhoa / (1 - rhoa - rhob)
  sigma2B <- rhob / (1 - rhoa - rhob)
  betaAB <-  gamma * sqrt(1 + sigma2A + sigma2B)
  ##  get predicted random effects
  if(get_effects)
    {
     a.est <- getEffects(rhoa,  list_eta = listeta1, list_w = listw1, niter = niter,
                         ws = ws, z = obj.gh$nodes) * sqrt(sigma2A)
     k1 <- length(unique(f1))
     a.out <- numeric(k1)
     a.out[sel1] <- 0
     a.out[setdiff(1:k1, sel1)] <- a.est
     b.est <- getEffects(rhob, list_eta = listeta2, list_w = listw2, niter = niter,
                         ws = ws, z = obj.gh$nodes) * sqrt(sigma2B)
     k2 <- length(unique(f2))
     b.out <- numeric(k2)
     b.out[sel2] <- 0
     b.out[setdiff(1:k2, sel2)] <- b.est
     }
  else a.out <- b.out <- NULL
  if(get_se)
  {
    df1f2 <- data.frame(f1 = f1, f2 = f2)
    vCL <- sandwich::vcovCL(obj_glm, cluster = df1f2, multi0 = TRUE)
    H <- (1 + sigma2A + sigma2B) * vCL
    beta.se <-  sqrt(diag(H))
  }
  else beta.se <- NULL
  return(list(beta = betaAB, beta.se = beta.se,
              sigmaA = sqrt(sigma2A), sigmaB = sqrt(sigma2B),
              a.est = a.out, b.est = b.out))
}


#' Fit an ordinal probit model with crossed random effects by the ARC method.
#' @param x Model matrix.
#' @param y Response vector (ordinal data with $K$ categories).
#' @param f1 Factor representing the row effects.
#' @param f2 Factor representing the columns effects.
#' @param obj_polr Object representing the fit of an ordinal probit model
#' as returned by \code{MASS::polr}.
#' @param nq Number of quadrature points. Default is 25.
#' @param niter Number of Newton-Raphson iteration to locate
#' the mode of each random effect. Default is 10.
#' @param get_effects Returns the random effects prediction.
#' Default is FALSE.
#' @param get_se Returns the standard errors of model coefficients.
#' Default is TRUE.
#' @return A list  containing the following components:
#' \item{\code{beta}}{the vector of estimated coefficients.}
#' \item{\code{beta.se}}{the estimated standard errors of the
#' coefficients.}
#' \item{\code{alpha}}{the vector of estimated thresholds.}
#' \item{\code{alpha.se}}{the estimated standard errors of the
#' thresholds.}
#' \item{\code{sigmaA, sigmaB}}{the estimated standard deviation of
#' row and column random effects.}
#' \item{\code{a.est, b.est}}{the predicted random effects.}
#' @export
#' @importFrom MASS polr
#' @import lme4
#' @import Rcpp
#' @examples
#' library(lme4)
#' mod.polr <- MASS::polr(factor(y) ~ studage + lectage + service + dept, method = "probit",
#'                        data = InstEval)
#' x <- model.matrix(mod.polr)[,-1]
#' fit.arc <- arcord.fit(x, InstEval$y, f1 = InstEval$s, f2 = InstEval$d, obj_polr = mod.polr)
#'
arcord.fit <- function(x, y, f1, f2, obj_polr, nq = 25, niter = 10,
                        get_effects = FALSE, get_se = TRUE){
  ## retrieve marginal parameters
  gamma <- obj_polr$coefficients
  zeta <- obj_polr$zeta
  alphae <- c(-100, zeta, 100)
  eta <- as.vector(x %*% gamma)

  ## now the row- and columns-steps
  listy1 <- split(y, f1)
  len1 <- unlist(lapply(listy1, length))
  sel1 <- which(len1 == 1)
  which.list1 <- which(len1 > 1)
  listy1 <- listy1[which.list1]
  listeta1 <- split(eta, f1)[which.list1]
  listy2 <- split(y, f2)
  len2 <- unlist(lapply(listy2, length))
  sel2 <- which(len2 == 1)
  which.list2 <- which(len2 > 1)
  listy2 <- listy2[which.list2]
  listeta2 <- split(eta, f2)[which.list2]

  ## obtain the tetrachoric correlations
  obj.gh <- statmod::gauss.quad(nq, "hermite")
  ws <- obj.gh$weights * exp(obj.gh$nodes^2)
  rhoa <- stats::optimise(likAGHOrd, interval = c(0, 1), alphae = alphae,
                   list_eta = listeta1, list_y = listy1, niter = niter,
                   ws = ws, z = obj.gh$nodes, maximum = TRUE)$maximum
  rhob <- stats::optimise(likAGHOrd, interval = c(0, 1), alphae = alphae,
                   list_eta = listeta2, list_y = listy2, niter = niter,
                   ws = ws, z = obj.gh$nodes, maximum = TRUE)$maximum
  ## back to the conditional parameterization
  sigma2A <- rhoa / (1 - rhoa - rhob)
  sigma2B <- rhob / (1 - rhoa - rhob)
  betaAB <-  gamma * sqrt(1 + sigma2A + sigma2B)
  alphaAB <- zeta * sqrt(1 + sigma2A + sigma2B)
  ##  get predicted random effects
  if(get_effects)
    {
     a.est <- getEffectsOrd(rhoa, alphae = alphae, list_eta = listeta1,
                            list_y = listy1, niter = niter, ws = ws,
                            z = obj.gh$nodes) * sqrt(sigma2A)
     k1 <- length(unique(f1))
     a.out <- numeric(k1)
     a.out[sel1] <- 0
     a.out[setdiff(1:k1, sel1)] <- a.est
     b.est <- getEffectsOrd(rhob, alphae = alphae, list_eta = listeta2,
                            list_y = listy2, niter = niter, ws = ws,
                            z = obj.gh$nodes) * sqrt(sigma2B)
     k2 <- length(unique(f2))
     b.out <- numeric(k2)
     b.out[sel2] <- 0
     b.out[setdiff(1:k2, sel2)] <- b.est
    }
    else a.est <- b.est <- NULL
  if(get_se)
  {
    df1f2 <- data.frame(f1 = f1, f2 = f2)
    vCL <- sandwich::vcovCL(obj_polr, cluster = df1f2, multi0 = TRUE)
    H <- (1 + sigma2A + sigma2B) * vCL
    se.all <- sqrt(diag(H))
    beta.se <-  se.all[1:length(gamma)]
    alpha.se <-  se.all[(length(gamma) + 1):length(se.all)]
  }
  else beta.se <- alpha.se <- NULL
  return(list(beta = betaAB, beta.se = beta.se,
              alpha = alphaAB, alpha.se = alpha.se,
              sigmaA = sqrt(sigma2A), sigmaB = sqrt(sigma2B),
              a.est = a.out, b.est = b.out))
}
