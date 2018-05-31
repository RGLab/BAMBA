##' Print a BAMBAResult Object
##'
##' This function prints basic information about the model fit by a
##' \code{\link{BAMBA}} call.
##'
##' @param x An object of class \code{BAMBAResult}.
##' @param ... Optional arguments; currently unused.
##' 
##' @export
##' 
##' @examples
##' print(bambaRes)
print.BAMBAResult <- function(x, ...) {
  nObs <- nrow(x$data)
  nSubj <- length(unique(x$data$subjectId))
  .catn <- function(n, name) {
      cat(n, ifelse(n == 1, name, paste0(name, "s")), "\n")
  }
  cat("A BAMBA model fit on",
      ifelse(x$dataType == "fc", "Fc array", "BAMA"), "data\n")
  .catn(nrow(x$data), "observation")
  .catn(length(unique(x$data$subjectId)), "subject")
  .catn(length(unique(x$data$groupId)), "group")
  .catn(length(unique(x$data$tp)), "timepoint")
  .catn(length(unique(x$data$agId)), "antigen")
  if (x$dataType == "fc") {
      .catn(length(unique(x$data$reId)), "reagent")
  }
}


##' Summarize a BAMBAResult Object
##'
##' This function prints basic information about the model fit by a
##' \code{\link{BAMBA}} call.
##'
##' @param object An object of class \code{BAMBAResult}.
##' @param ... Optional arguments; currently unused.
##' 
##' @export
##' 
##' @examples
##' summary(bambaRes)
summary.BAMBAResult <- function(object, ...) {
    print(object, ...)
}
