#' Gets the parameter summary from a stanfit object
#'
#' @param modelFit The stanfit object to summarize.
#'
#' @return A \code{data.frame} containing summaries of all
#'   sampled parameters in the model.
#'
#' @noRd
build_BAMBA_summary <- function(modelFit) {
    parSumm <- as.data.frame(rstan::summary(modelFit)$summary)
    parSumm$var <- rownames(parSumm)
    names(parSumm) <- c("mean", "se_mean", "sd",
                         "q025", "q25", "q50", "q75", "q975",
                       "n_eff", "Rhat", "var")

    parSumm %>% dplyr::select(var, mean, se_mean, sd, q025, q25, q50, q75, q975, n_eff, Rhat)
}
