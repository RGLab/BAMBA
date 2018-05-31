#' Calculate the expected false-discovery rate for a given response cutoff
#'
#' This function calculates the expected false-discovery rate for a
#' set of response probabilities given the probability cutoff for
#' calling a response.
#'
#' @param responseProbs The response probabilities to use for calculating
#'   the expected false-discovery rate.
#' @param cutoff The probability cutoff probability for calling a response.
#' 
#' @return The expected false-discovery rate.
#'
#' @export
efdr <- function(responseProbs, cutoff = 0.5) {
    probs <- responseProbs[responseProbs >= cutoff]
    sum(1-probs) / length(probs)
}
