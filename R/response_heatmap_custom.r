#' Plot the response probabilities for all observations as a heatmap
#'
#' This function plots the response probabilities from a model fit as a
#' heatmap with additional options for sorting and filtering both axes.
#' 
#' @param result The BAMBAResult object.
#' @param xorderTable A \code{data.frame} with all ag/re/tp combinations
#'   to include as well as ordering, labeling, and color information.
#'   Should have the following columns: ag, re, tp, order, label, color.
#' @param yorderTable A \code{data.frame} with all subjectIds to include
#'   as well as ordering, labeling, and color information.
#'   Should have the following columns: subjectId, order, label, color
#' @param responseThreshold If not NULL, the threshold probability
#'   defining a response, resulting in a two-color heatmap rather than
#'   a continuous heatmap. Defaults to \code{NULL}
#' @param xtext The label for the x-axis. Defaults to 'Antigen/Fc Variable'.
#' @param xlines A string defining the color for lines separating groups
#'   (by label) on the x-axis or \code{NULL} for no lines.
#'   Defaults to 'white'.
#' @param ytext The label for the y-axis. Defaults to 'SubjectId'.
#' @param ylines A string defining the color for lines separating groups
#'   (by label) on the y-axis or \code{NULL} for no lines.
#'   Defaults to \code{NULL}
#'
#' @return A ggplot heatmap.
#'
#' @export
response_heatmap_custom <- function(result,
                             xorderTable,
                             yorderTable,
                             responseThreshold = NULL,
                             xtext = "Antigen/Fc Variable",
                             xlines = "white",
                             ytext = "SubjectId",
                             ylines = NULL) {

    resp <- responses(result)
    if (!"re" %in% names(resp)) {
        resp$re <- "noFc"
        xorderTable$re <- "noFc"
    }

    resp <- resp %>%
        dplyr::filter(paste(ag, re, tp) %in%
                      paste(xorderTable$ag, xorderTable$re, xorderTable$tp),
                      subjectId %in% yorderTable$subjectId)
    hmData <- resp %>%
        left_join(xorderTable %>%
                  rename(xorder = order,
                         xlabel = label,
                         xgroup = group,
                         xcolor = color),
                  by = c("ag", "re", "tp")) %>%
        left_join(yorderTable %>%
                  rename(yorder = order,
                         ylabel = label,
                         ygroup = group,
                         ycolor = color),
                  by = "subjectId") %>%
        mutate(xorder = dense_rank(xorder),
               yorder = dense_rank(yorder)) %>%
        arrange(xorder, yorder)

    xax <- hmData %>%
        dplyr::select(xorder, xlabel, xgroup, xcolor) %>%
        distinct() %>%
        arrange(xorder)

    xsep <- cumsum(rle(xax$xlabel)$lengths)
    ##    xsep2 <- floor((xsep + c(0, xsep[1:(length(xsep)-1)])) / 2)
    xsep2 <- NULL
    if (length(xsep) > 1) {
        xsep2 <- ceiling((xsep + c(0, xsep[1:(length(xsep)-1)])) / 2)
    } else {
        xsep2 <- ceiling(xsep /2)
    }
    ## xgrpsep <- cumsum(table(xax$xgroup))
    xgrpsep <- cumsum(rle(xax$xgroup)$lengths)


    yax <- hmData %>%
        dplyr::select(yorder, ylabel, ygroup, ycolor) %>%
        distinct() %>%
        arrange(yorder)
    ysep <- cumsum(rle(yax$ylabel)$lengths)
##    ysep2 <- floor((ysep + c(0, ysep[1:(length(ysep)-1)])) / 2)
    ysep2 <- NULL
    if (length(ysep) > 1) {
        ysep2 <- ceiling((ysep + c(0, ysep[1:(length(ysep)-1)])) / 2)
    } else {
        ysep2 <- ceiling(ysep /2)
    }
    ## ygrpsep <- cumsum(table(yax$ygroup))
    ygrpsep <- cumsum(rle(yax$ygroup)$lengths)

    hmPlot <- NULL
    if (is.null(responseThreshold)) {
        hmPlot <- ggplot(hmData) +
            geom_tile(aes(xorder, yorder, fill=responseProb))
    } else {
        hmPlot <- ggplot(hmData) +
            geom_tile(aes(xorder, yorder,
                          fill=responseProb > responseThreshold)) +
            scale_fill_manual(values = c("darkred", "darkgreen"))
    }
    
    if (!is.null(xlines)) {
        hmPlot <- hmPlot +
            geom_vline(xintercept = xgrpsep+0.5, color=xlines)
    }
    if (!is.null(ylines)) {
        hmPlot <- hmPlot +
            geom_hline(yintercept = ygrpsep+0.5, color=ylines)
    }

    hmPlot +
        scale_x_discrete(xtext, limits=xax$xorder,
                         breaks=xax$xorder[xsep2], labels=xax$xlabel[xsep2]) +
        scale_y_discrete(ytext, limits=yax$yorder,
                         breaks=yax$yorder[ysep2], labels=yax$ylabel[ysep2]) +
        theme(axis.text.x = element_text(angle = 45, hjust=.5, vjust=.5,
                                         color = xax$xcolor[xsep2]),
              axis.text.y = element_text(color = yax$ycolor[ysep2]))
}
