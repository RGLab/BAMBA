#' @include utils.r
#' @include data_maps.r

#' Check and prepares data for use in BAMBA modeling
#'
#' This function validates and prepares data for use by stan in
#' fitting the BAMBA model.
#'
#' @param data The data to validate and prepare.
#'
#' @return A validated and prepared \code{data.frame}, ready for
#'   use in the BAMBA model.
#'
#' @noRd
prepare_data <- function(data) {
    errorStrings <- NULL

    ## Check for the appropriate columns in the data
    colErrorStr <- paste0("Invalid data - must contain the following columns: ",
                         paste(BAMBA_COLUMN_NAMES, collapse = ", "), "\n")
    if ( !all( BAMBA_COLUMN_NAMES %in% colnames(data) ) ) {
        errorStrings <- c(errorStrings, colErrorStr)
    }

    ## Ensure group, re, tp columns for fc array data
    ## Add them for bama data as well, but silently
    if (!"group" %in% colnames(data)) {
      message("data - simple group column added")
      data <- data %>% mutate(group = "No Group")
    }
    if (!"re" %in% colnames(data)) {
      message("data - simple re column added")
      data <- data %>% mutate(re = "No Fc Variable")
    }
    if (!"tp" %in% colnames(data)) {
      message("data - simple tp column added")
      data <- data %>% mutate(tp = 1)
    }
    if (!"sat" %in% colnames(data)) {
      message("data - simple sat column added")
      data <- data %>% mutate(sat = FALSE)
    }
    
    ## Check tp column
    data$tp <- suppressWarnings(as.numeric(data$tp))
    if (any(is.na(data$tp))) {
        errorStrings <- c(errorStrings,
                          "Invalid data - tp column must be numeric\n")
    }

    ## Check sat column
    data$sat <- suppressWarnings(as.logical(data$sat))
    if (any(is.na(data$sat))) {
      errorStrings <- c(errorStrings,
                        "Invalid data - sat column must contain TRUE or FALSE values only\n")
    }
    
    if(!is.null(errorStrings)) {
        stop(errorStrings)
    }

    ## Standardize tp column
    tpRank <- dense_rank(data$tp)
    if (any(data$tp != tpRank)) {
        warning("data - tp column transformed to dense ranks")
        data$tp <- tpRank
    }

    ## Check val column
    data$val <- suppressWarnings(as.numeric(data$val))
    nNAvals = sum(is.na(data$val))
    if (nNAvals > 0) {
        warning(str_c(nNAvals, " non-numeric or NA value",
                      ifelse(nNAvals == 1, "", "s"),
                      " found - removing for analysis"))
    }

    ## Return prepared data
    data %>%
        dplyr::filter(!is.na(val)) %>%
        mutate(subjectId = as.character(subjectId),
               group = as.character(group),
               ag = as.character(ag),
               re = as.character(re)) %>%
        left_join(group_map(data), by = "group") %>%
        left_join(antigen_map(data), by = "ag") %>%
        left_join(reagent_map(data), by = "re") %>%
        arrange(group, subjectId, ag, re, tp)

}
