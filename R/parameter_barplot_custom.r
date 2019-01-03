#' Plot the parameter data as a barplot
#'
#' This function plots sampled parameter values from a model fit as a
#' barplot. All sorting and filtering is done in the passed data.
#'
#' @param paramData A \code{data.frame} with all the data necessary for
#' creating the barplot. This includes the following columns:
#' \describe{
#' \item{order}{The order of parameters on the axis. Should include all
#' integers from 1 to the number of rows in the dataset.}
#' \item{label}{The parameter label to use in the barplot. Consecutive
#' parameters with the same label will share a label on the axis rather
#' than repeating the label.}
#' \item{labelColor}{The color to use for the label text.}
#' \item{paramValue}{The value of the parameter (length of the bar).}
#' \item{qlower}{The lower sampling percentile for the value of the
#' parameter. Required if \code{incIntervals = TRUE}.}
#' \item{qupper}{The upper sampling percentile for the value of the
#' parameter. Required if \code{incIntervals = TRUE}.}
#' \item{group}{The parameter group. Each distinct group should cover
#' a contiguous set of parameters by order in order to separate groups
#' by lines. Required if \code{groupSepColor != NULL}.}
#' \item{fillGroup}{The group designation for bar fill color.}
#' \item{fillColor}{The fill color for the group.}
#' }
#' @param paramLabel The label for the parameter axis. Defaults to 'Parameter'.
#' @param valueLabel The label for the value axis. Defaults to 'Value'.
#' @param fillLabel The label for the fill color legend. Defaults to 'Fill'.
#' @param groupSepColor A string defining the color for lines separating groups
#' on the parameter axis or \code{NULL} for no lines. Defaults to \code{NULL}.
#' @param incIntervals A boolean indicating whether to include sampling
#' intervals for the parameter values. Defaults to \code{FALSE}.
#' @param intervalPtSize The size for the center point of the
#' sampling interval. Defaults to 3.
#' @param intervalBarSize The size for the sampling interval bars.
#' Defaults to 1.5
#' @param intervalColor A string defining the color for the sampling interval
#' error bars. Defaults to 'gray30'.
#'
#' @return A ggplot barplot object.
#'
#' @export
parameter_barplot_custom <- function(paramData,
                                     paramLabel = "Parameter",
                                     valueLabel = "Value",
                                     fillLabel = "Fill",
                                     groupSepColor = NULL,
                                     incIntervals = F,
                                     intervalPtSize = 3,
                                     intervalBarSize = 1.5,
                                     intervalColor = "gray30") {
    paramData <- paramData %>%
        mutate(fillGroup = as.character(fillGroup))
    barColors = paramData %>%
        arrange(order) %>%
        dplyr::select(fillGroup, fillColor) %>%
        distinct()
    xax <- paramData %>%
        dplyr::select(order, label, group, labelColor) %>%
        distinct() %>%
        arrange(order)
    xsep <- cumsum(table(xax$label))
    xsep2 <- ceiling((xsep + c(0, xsep[1:(length(xsep)-1)])) / 2)
    xgrpsep <- cumsum(table(xax$group))
        
    barPlot <- ggplot(paramData, aes(order, paramValue)) +
        labs(y = valueLabel) +
        coord_flip()
    if (nrow(barColors) > 1) {
        barPlot <- barPlot +
            geom_bar(aes(fill = fillGroup),
                     stat = "identity", alpha = 0.9) +
        scale_fill_manual(limits = barColors$fillGroup,
                          values = barColors$fillColor) +
            labs(fill = fillLabel)
    } else {
        barPlot <- barPlot +
            geom_bar(stat = "identity", alpha = 0.9)
    }
    if (!is.null(groupSepColor)) {
        barPlot <- barPlot +
            geom_vline(xintercept = xgrpsep+0.5, color=groupSepColor)
    }
    if (incIntervals) {
        barPlot <- barPlot +
            geom_point(size = intervalPtSize,
                       color = intervalColor) +
            geom_errorbar(aes(ymin = qlower,
                              ymax = qupper),
                          size = intervalBarSize,
                          color = intervalColor)
    }
    barPlot +
        scale_x_discrete(paramLabel, limits = xax$order,
                         breaks = xax$order[xsep2],
                         labels = xax$label[xsep2]) +
        theme(axis.text.x = element_text(color = xax$labelColor[xsep2]))
}
