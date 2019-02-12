#' Fit the BAMBA model on BAMA or Fc array data
#'
#' This function fits the \code{BAMBA} model.
#' TODO: change params out to just list all the options?
#'       saving the stanfit object or returning it
#'       keeping original data
#'
#' @param data The data to be modeled.
#' @param nChains The number of chains to run for the stan sampling.
#'   Defaults to 1. If more than one chain is used, it is recommended
#'   to call \code{options(mc.cores = parallel::detectCores())}
#'   to allow for parallel processing of multiple chains.
#' @param nIter The number of iterations per chain for the stan sampling.
#'   Defaults to 2000.
#' @param outFolder The folder to save the stan results to.
#'   If disk space is an issue, you may set this to \code{NULL} to
#'   not save the stanfit object. Defaults to \code{NULL}
#' @param outFile The filename to save the stan results to.
#'   Defaults to '<yyyy-mm-dd>_BAMBA_stanfit.rds'.
#' @param ... Additional parameters to pass to the stan sampling function.
#'   See the rstan documentation for more information.
#'
#' @return A \code{BAMBAResult} is a list with the following components:
#'
#' \item{data}{The validated and prepared data used to fit
#'   the \code{BAMBA} model.}
#' \item{chains}{The number of chains run for the stan sampling.}
#' \item{iter}{The number of iterations per chain for the stan sampling.}
#' \item{parameters}{A \code{data.frame} summarizing all parameters sampled
#'   for the model.}
#' \item{mu0}{A \code{data.frame} summarizing the samples of the
#'   baseline mean parameter, mu0.}
#' \item{mu_ag}{A \code{data.frame} summarizing the samples of the
#'   antigen offsets, mu_ag.}
#' \item{mu_re}{A \code{data.frame} summarizing the samples of the
#'   Fc variable offsets, mu_re. Only included if \code{dataType == 'fc'}.}
#' \item{mu_ar}{A \code{data.frame} summarizing the samples of the
#'   antigen/Fc variable offsets, mu_ar = mu_ag + mu_re.
#'   mu_ar has more sampling stability than mu_ag and mu_re.}
#' \item{omega_t}{A \code{data.frame} summarizing the samples of the
#'   prior response probabilities per timepoint, omega_t.}
#' \item{omega_ag}{A \code{data.frame} summarizing the samples of the
#'   prior response probabilities per antigen, omega_ag.}
#' \item{omega_re}{A \code{data.frame} summarizing the samples of the
#'   prior response probabilities per Fc variable, omega_re.}
#' \item{omega_grp}{A \code{data.frame} summarizing the samples of the
#'   prior response probabilities per group, omega_grp.}
#' \item{hyperparameters}{A \code{data.frame} summarizing
#'   the samples of the model hyperparameters.}
#'
#' @export
#' @example examples/BAMBA_fit.r
BAMBA <- function(data,
                  nChains = 1,
                  nIter = 2000,
                  outFolder = NULL,
                  outFile = date_filename("BAMBA_stanfit.rds"),
                  ...) {

    ## Model preparation
    data <- prepare_data(data)
    modelData <- build_model_data(data)
    paramInit <- build_parameter_initialization(modelData)
    modelName <- build_model_name(modelData)

    ## Run stan
    stanRes <- rstan::sampling(stanmodels[[modelName]],
                    data = modelData,
                    init = function(){paramInit},
                    iter = nIter,
                    chains = nChains,
                    include = FALSE,
                    pars = c("ystar", "soft_z"),
                    ...)

    if (!is.null(outFolder)) {
        saveRDS(stanRes, file.path(outFolder, outFile))
    }
    
    output <- list(
        data = data,
        chains = nChains,
        iter = nIter,
        parameters = build_BAMBA_summary(stanRes))

    output <- add_parameter_summaries(output)

    class(output) <- "BAMBAResult"

    return(output)
}
