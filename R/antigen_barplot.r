#' Plot the antigen means as a barplot
#'
#' This function plots the antigen means from a model fit as a
#'   barplot.
#'
#' @param result The BAMBAResult object.
#' @param includeIntervals A boolean indicating whether to
#'   include sampling intervals for the antigen means.
#'   Defaults to \code{FALSE}.
#'
#' @return A ggplot barplot object.
#'
#' @export
antigen_barplot <- function(result,
                            includeIntervals = F) {
    agData <- result$mu_ag %>%
        dplyr::select(ag, tp, mu_ag, q025, q975) %>%
        arrange(ag, tp) %>%
        rename(paramValue = mu_ag,
               qlower = q025,
               qupper = q975) %>%
        mutate(order = 1:n(),
               label = ag,
               group = ag,
               labelColor = "black",
               fillGroup = tp,
               fillColor = rainbow(n=length(unique(tp)),
                                   s = 0.8)[tp])

    parameter_barplot_custom(agData,
                             paramLabel   = "Antigen",
                             valueLabel   = "Mean Offset",
                             fillLabel    = "Timepoint",
                             incIntervals = includeIntervals)
}
