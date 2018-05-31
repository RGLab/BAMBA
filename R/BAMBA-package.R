#' BAMBA (Bayesian Analysis of Multiplexed Binding Antibody data)
#'
#' This package fits a Bayesian grouped mixture model for BAMA and
#' Fc array data using stan. The model identifies responses in the
#' data.
#'
#' @docType package
#' @name BAMBA-package
#' @aliases BAMBA-package
#' @useDynLib BAMBA, .registration = TRUE
#'
#' @import methods
#' @import stats
#' @import Rcpp
#' @import rstantools
#' @import stringr
#' @import tidyr
#' @import ggplot2
#'
#' @seealso
#'   \itemize{
#'     \item \code{\link{BAMBA}}, for the main model fitting routine.
#'   }
#'
#' 
NULL
