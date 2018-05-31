#' Calculate the cutoff probability for a requested false-discovery rate
#'
#' This function calculates the minimum cutoff probability for defining a
#' response from a set of response probabilities such that the expected
#' false-discovery rate is no greater than the requested FDR.
#'
#' @param responseProbs The response probabilities to use for calculating
#'   the minimum cutoff probability.
#' @param fdr The requested false-discovery rate.
#' 
#' @return The maximum probability cutoff for defining a response that
#'   will result in an expected false-discovery rate equal to or less
#'   than that requested.
#'
#' @export
fdr_cutoff <- function(responseProbs, fdr = 0.05) {
    probs <- sort(responseProbs, decreasing = TRUE)
    efdr <- cumsum(1-probs) / (1:length(probs))
    min(probs[efdr <= fdr])
}
