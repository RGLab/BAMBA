#' Get the inferred response probabilities for all observations from
#' a model fit
#'
#' This function adds response probabilities from a model fit
#' to the modeled dataset.
#'
#' @param result The BAMBAResult object to query.
#'
#' @return A \code{data.frame} containing the observation and response
#'   probability information.
#'
#' @export
responses <- function(result) {
    respProbs <- result %>% extract_params2("z[", ",2]") %>%
        rename(responseProb = mean) %>%
        mutate(isResp = responseProb > 0.5) %>%
        select(isResp, responseProb, sd, q025, q975, n_eff, Rhat)
    ## make sure to grab the correct columns
    bind_cols(result$data, respProbs) %>%
        select(subjectId, group, groupId, ag, agId,
               re, reId, tp, val, isResp, responseProb,
               sd, q025, q975, n_eff, Rhat)
}
