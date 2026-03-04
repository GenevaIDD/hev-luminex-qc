# R/utils.R
library(stringr)

guess_plate_id <- function(path) {
  # Example: PlateRunResults_PLATE_02032026_plate6.csv
  bn <- basename(path)
  m <- str_match(bn, "PlateRunResults_([^\\.]+)\\.csv")[,2]
  ifelse(is.na(m), tools::file_path_sans_ext(bn), m)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))