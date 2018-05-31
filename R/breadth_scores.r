.breadth_scores <- function(resp) {
    fab <- resp %>%
        group_by(group, subjectId, ag) %>%
        summarize(prob = max(responseProb)) %>%
        group_by(group, subjectId) %>%
        summarize(fabBreadth = mean(prob))
    fcr <- resp %>%
        group_by(group, subjectId, re) %>%
        summarize(prob = max(responseProb)) %>%
        group_by(group, subjectId) %>%
        summarize(fcRBreadth = mean(prob))
    fab %>%
        left_join(fcr, by = c("group", "subjectId"))
}

#' Gets fab and FcR breadth scores per subject from
#' a BAMBAResult object
#'
#' @param result The BAMBAResult object.
#' @param agClasses A named list of antigen classes, each item being
#' a vector of antigen names used in the data to filter by before computing scores.
#' Defaults to \code{NULL}, indicating no filtering.
#' @param reClasses A named list of Fc variable classes, each item being
#' a vector of Fc variable names used in the data to filter by before computing scores.
#' Defaults to \code{NULL}, indicating no filtering.
#' @param includeMagBreadth A boolean indicating whether to include
#' magnitude-weighted breadth scores. Defaults to \code{FALSE}.
#'
#' @return A \code{data.frame} containing breadth scores for each subject,
#' for all antigens as well as for each antigen class included in agClasses.
#'
#' @export
breadth_scores <- function(result,
                           tps = unique(result$data$tp),
                           agClasses = NULL,
                           reClasses = NULL) {
    scores <- .collect_scores(result,
                              tps,
                              agClasses,
                              reClasses,
                              .breadth_scores)
    if (result$dataType == "bama") {
        scores %>% select(-reClass, -fcRBreadth)
    }
    else { scores }
}

