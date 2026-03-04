# R/qc_standards.R
library(dplyr)
library(stringr)
library(tidyr)

extract_standards <- function(df) {
  df %>%
    filter(sample_type %in% c("standard", "blank")) %>%
    mutate(
      # Pull "1/200" from sample like "Std 1/200"
      std_dilution = str_match(sample, "(\\d+\\s*/\\s*\\d+)")[,2] %>%
        str_replace_all("\\s+", ""),
      denom = suppressWarnings(as.numeric(str_match(std_dilution, "/\\s*(\\d+)")[,2])),
      is_blank = sample_type == "blank"
    ) %>%
    # Keep only true standards (blanks have denom NA)
    mutate(
      x = if_else(!is_blank & !is.na(denom), log10(denom), NA_real_)
    )
}

standard_summary <- function(std_df) {
  std_df %>%
    filter(!is_blank, !is.na(denom)) %>%
    group_by(plate_id, analyte) %>%
    summarise(
      n_levels = n_distinct(denom),
      min_denom = min(denom, na.rm = TRUE),
      max_denom = max(denom, na.rm = TRUE),
      .groups = "drop"
    )
}