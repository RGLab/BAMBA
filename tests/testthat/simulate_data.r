library(dplyr)
library(stringr)

#' Simulate data for Fc modeling
#'
#' Simulates Fc array data according to given parameters
#' for testing Fc modeling. Allows the variance to differ
#' between ag/re combinations
#'
#' @param nTp number of timepoints
#' @param nGrp number of groups
#' @param nSubj number of subjects per group
#' @param nAg number of antigens in the array
#' @param nRe number of reagents in the array
#' @param agAlpha alpha parameter for randomly generating
#' antigen response means from a beta distribution
#' @param agBeta beta parameter for randomly generating
#' antigen response means from a gamma distribution
#' @param reAlpha alpha parameter for randomly generating
#' reagent response means from a gamma distribution
#' @param reBeta beta parameter for randomly generating
#' reagent response means from a gamma distribution
#' @param sigmaAlpha alpha parameter for randomly generating
#' ag/re group variances from a gamma distribution
#' @param sigmaBeta beta parameter for randomly generating
#' ag/re group variances from a gamma distribution
#' @param mu0 non-reponse mean
#' @param groupProbs prior group response probabilities
#' @param agProbs prior antigen response probabilities
#' @param reProbs prior reagent response probabilities
#'
#' @return a list containing the simulated data as well as
#' the antigen and reagent response meanse and the ag/re group variances
#'
#' @export
simulate_fc_data <- function(nTp,
                             nGrp, nSubj,
                             nAg, nRe,
                             agAlpha, agBeta,
                             reAlpha, reBeta,
                             sigmaAlpha, sigmaBeta,
                             mu0, groupProbs,
                             agProbs, reProbs) {

    N <- nTp * nGrp * nSubj * nAg * nRe
    muAg <- rgamma(nAg, agAlpha, agBeta)
    muRe <- rgamma(nRe, reAlpha, reBeta)
    sigma <- rgamma(nAg * nRe, sigmaAlpha, sigmaBeta)
    simData <- expand.grid(1:nGrp, 1:nSubj, 1:nAg, 1:nRe, 1:nTp) %>%
        rename(group = Var1,
               subjectId = Var2,
               ag = Var3,
               re = Var4,
               tp = Var5) %>%
        mutate(responseProb = groupProbs[group] * agProbs[ag] * reProbs[re],
               isResponse = rbinom(N, 1, prob = responseProb),
               val = rnorm(N) * sigma[ag + (re-1)*nAg] + mu0 + isResponse * (muAg[ag] + muRe[re]),
               group = paste0("grp", str_pad(group, 3, "left", "0")),
               subjectId = paste0(group, "_subj", str_pad(subjectId, 3, "left", "0")),
               ag = paste0("ag", str_pad(ag, 3, "left", "0")),
               re = paste0("re", str_pad(re, 3, "left", "0")))

    list(simData = simData,
         agMeans = muAg,
         reMeans = muRe,
         sigmas = sigma)
}
