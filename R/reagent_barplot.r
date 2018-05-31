#' Plot the reagent means as a barplot
#'
#' This function plots the reagent means from a model fit as a
#'   barplot.
#'
#' @param result The BAMBAResult object.
#' @param includeIntervals A boolean indicating whether to
#'   include sampling intervals for the reagent means.
#'   Defaults to \code{FALSE}.
#'
#' @return A ggplot barplot object.
#'
#' @export
reagent_barplot <- function(result,
                            includeIntervals = F) {
    reData <- result$mu_re %>%
        select(re, tp, mu_re, q025, q975) %>%
        arrange(re, tp) %>%
        rename(paramValue = mu_re,
               qlower = q025,
               qupper = q975) %>%
        mutate(order = 1:n(),
               label = re,
               group = re,
               labelColor = "black",
               fillGroup = tp,
               fillColor = rainbow(n=length(unique(tp)),
                                   s = 0.8)[tp])

    parameter_barplot_custom(reData,
                             paramLabel   = "Fc Variable",
                             valueLabel   = "Mean Offset",
                             fillLabel    = "Timepoint",
                             incIntervals = includeIntervals)
}
