# R/qc_history.R
library(dplyr)

history_load <- function(path) {
  if (file.exists(path)) readRDS(path) else tibble()
}

history_append <- function(history_df, new_df, key_cols) {
  # Avoid duplicates by key
  out <- bind_rows(history_df, new_df) %>%
    distinct(across(all_of(key_cols)), .keep_all = TRUE)
  out
}

history_save <- function(df, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  saveRDS(df, path)
}

make_history_records <- function(df, rep_tbl, blank_tbl) {
  # Per plate/analyte replicate fail rate
  rep_hist <- rep_tbl %>%
    mutate(fail = flag == "fail_high_cv") %>%
    summarise(
      n_pairs = n(),
      pct_fail = mean(fail, na.rm = TRUE),
      .by = c(plate_id, run_date, analyte)
    ) %>%
    mutate(metric = "replicate_pct_fail")
  
  blank_hist <- blank_tbl %>%
    transmute(
      plate_id,
      analyte,
      metric = "blank_q99",
      value = blank_q99
    )
  
  list(rep_hist = rep_hist, blank_hist = blank_hist)
}