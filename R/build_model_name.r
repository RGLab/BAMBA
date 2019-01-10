build_model_name <- function(modelData) {
  mName <- "bamba_ag"
  if (modelData$N_re > 1) {
    mName <- paste(mName, "fc", sep="_")
  }
  if (modelData$N_grp > 1) {
    mName <- paste(mName, "grp", sep="_")
  }
  if (modelData$N_sat > 0) {
    mName <- paste(mName, "sat", sep="_")
  }
  return(paste(mName, "model", sep="_"))
}