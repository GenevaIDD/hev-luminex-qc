# R/qc_standards_history.R
library(dplyr)
library(stringr)

extract_std_points <- function(df) {
  # Standards only; keep replicates separately
  df %>%
    filter(sample_type == "standard") %>%
    mutate(
      std_dilution = str_match(sample, "(\\d+\\s*/\\s*\\d+)")[,2] %>% str_replace_all("\\s+", ""),
      denom = suppressWarnings(as.numeric(str_match(std_dilution, "/\\s*(\\d+)")[,2]))
    ) %>%
    filter(!is.na(denom)) %>%
    select(plate_id, run_date, analyte, denom, value, rep_index) %>%
    mutate(run_date = as.character(run_date))
}

std_history_load <- function(path) {
  if (file.exists(path)) readRDS(path) else tibble::tibble()
}

std_history_save <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
}

std_history_append <- function(hist, new_pts) {
  # keep unique by plate/analyte/denom/rep_index
  dplyr::bind_rows(hist, new_pts) %>%
    distinct(plate_id, analyte, denom, rep_index, .keep_all = TRUE)
}