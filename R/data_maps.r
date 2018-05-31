.uniqsort <- function(x) {
    x %>% unique() %>% sort() %>% as.character()
}

group_map <- function(data) {
    grps <- data$group %>% .uniqsort()
    data_frame(group = grps,
               groupId = 1:length(grps))
}

antigen_map <- function(data) {
    ags <- data$ag %>% .uniqsort()
    data_frame(ag = ags,
               agId = 1:length(ags))
}

reagent_map <- function(data) {
    res <- data$re %>% .uniqsort()
    data_frame(re = res,
               reId = 1:length(res))
}
