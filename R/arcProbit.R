#' arcProbit: All-row-col fitting of categorical data models with crossed random effects
#'
#'
#' The \code{arcProbit} package can be used to estimate the parameters of mixed models for binary
#' and ordinal data with crossed random effects. It implements the ARC estimation proposed in
#' https://arxiv.org/abs/2308.15681. A vignette is available with some illustrative examples.
#'
#'
#' @docType package
#' @name arcProbit
#' @references Bellio, R., Ghosh, S., Owen, A. and Varin, C. (2023). Scalable Estimation of Probit Models with
#' Crossed Random Effects. https://arxiv.org/abs/2308.15681
#' @useDynLib binCross
#' @useDynLib ordCross
#' @useDynLib arcProbit, .registration = TRUE
NULL
#> NULL
