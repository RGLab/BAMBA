.collect_scores <- function(result,
                            tps,
                            agClasses,
                            reClasses,
                            breadthFn,
                            ...) {
    resp <- responses(result) %>%
        dplyr::filter(tp %in% tps)
    tpLabel = paste(tps, collapse = ", ")
    scores <- breadthFn(resp, ...) %>%
        mutate(tp = tpLabel,
               agClass = "--All--",
               reClass = "--All--")
    for (r in names(reClasses)) {
        rScores <- breadthFn(
            resp %>% dplyr::filter(re %in% reClasses[[r]]),
            ...
        ) %>% mutate(tp = tpLabel,
                     agClass = "--All--",
                     reClass = r)
        scores <- scores %>%
            bind_rows(rScores)
    }
    for (c in names(agClasses)) {
        cScores <- breadthFn(
            resp %>% dplyr::filter(ag %in% agClasses[[c]]),
            ...
        ) %>% mutate(tp = tpLabel,
                     agClass = c,
                     reClass = "--All--")
        scores <- scores %>%
            bind_rows(cScores)
        for (r in names(reClasses)) {
            rScores <- breadthFn(
                resp %>% dplyr::filter(ag %in% agClasses[[c]],
                                       re %in% reClasses[[r]]),
                ...
            ) %>% mutate(tp = tpLabel,
                         agClass = c,
                         reClass = r)
            scores <- scores %>%
                bind_rows(rScores)
        }
    }
    if (length(tps) > 1) {
        for (t in tps) {
            scores <- scores %>%
                bind_rows(.collect_scores(result, t, agClasses,
                                          reClasses, breadthFn, ...))
        }
    }
    scores
}
