date_filename <- function(filename) {
    paste(format(Sys.Date(), "%Y-%m-%d"),
          filename, sep="_")
    
}
