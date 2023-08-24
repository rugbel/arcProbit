#' arcProbit: All-row-col fitting of categorical data models with crossed random effects
#'
#'
#' The \code{arcProbit} package can be used to estimate the parameters of a p2 model
#' for directed binary networks with correlated random effects. It implements
#' (approximate) maximum likelihood estimation, following the methodology
#' studied in ARC REF.
#'
#' A vignette is available with some illustrative examples.
#'
#'
#' @docType package
#' @name arcProbit
#' @references Bellio, R., Ghosh, S., Owen, A. and Varin, C. (2023). Scalable Estimation of Probit Models with
#' Crossed Random Effects. Manuscript in preparation.
#' @useDynLib binCross
#' @useDynLib ordCross
#' @useDynLib arcProbit, .registration = TRUE
NULL
#> NULL
