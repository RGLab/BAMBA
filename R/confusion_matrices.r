##' Calculate a confusion matrix from response data
##'
##' This function calculates a vector of cutoffProb = cutoff,
##' tp, fp, tn, and fn according to provided response data.
##'
##' @param cutoff The cutoff probability to use for defining a response.
##' @param rData The response data to check. Must include columns
##' isResponder and responseProb.
##'
##' @return A vector containing the given cutoff probability and
##' confusion matrix entries tp, fp, tn, and fn.
##'
##' @noRd
.conf <- function(cutoff, rData) {
    c(cutoff, sum(rData$responseProb >= cutoff & rData$isResponder),
      sum(rData$responseProb >= cutoff & !rData$isResponder),
      sum(rData$responseProb < cutoff & !rData$isResponder),
      sum(rData$responseProb < cutoff & rData$isResponder))
}

##' Generate a table of confusion matrices response data
##'
##' @param respData The response data to check. Must include columns
##' isResponder and responseProb.
##' @param cutoffs A set of cutoff probabilities for defining response.
##'
##' @return A \code{data.frame} of confusion matrices for the provided
##' cutoff probabilities with columns cutoffProb, tp, fp, tn, and fn.
##'
##' @noRd
.confusion_matrices <- function(respData, cutoffs) {
    confMat <- as.data.frame(t(sapply(cutoffs, .conf, respData)))
    names(confMat) <- c("cutoffProb", "tp", "fp", "tn", "fn")
    confMat
}
