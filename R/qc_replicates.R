# R/qc_replicates.R
library(dplyr)
library(tidyr)

compute_replicate_variability <- function(df,
                                          cv_threshold = 0.30,
                                          log_pseudocount = 1,
                                          blank_quantile = 0.99,
                                          low_signal_rule = c("blank_q99", "fixed"),
                                          low_signal_fixed = 50) {
  
  low_signal_rule <- match.arg(low_signal_rule)
  
  # blank cutoff per plate/analyte
  blank_cutoffs <- df %>%
    filter(sample_type == "blank") %>%
    group_by(plate_id, analyte) %>%
    summarise(blank_q = quantile(value, blank_quantile, na.rm = TRUE), .groups = "drop")
  
  reps <- df %>%
    filter(sample_type == "specimen") %>%
    group_by(plate_id, run_date, plate_batch, sample, dilution, analyte) %>%
    summarise(
      n = sum(!is.na(value)),
      rep1 = ifelse(n >= 1, value[which(!is.na(value))[1]], NA_real_),
      rep2 = ifelse(n >= 2, value[which(!is.na(value))[2]], NA_real_),
      mean_raw = mean(value, na.rm = TRUE),
      sd_raw   = sd(value, na.rm = TRUE),
      cv_raw   = sd_raw / mean_raw,
      mean_log = mean(log10(value + log_pseudocount), na.rm = TRUE),
      sd_log   = sd(log10(value + log_pseudocount), na.rm = TRUE),
      cv_log   = sd_log / mean_log,
      .groups = "drop"
    ) %>%
    left_join(blank_cutoffs, by = c("plate_id", "analyte")) %>%
    mutate(
      low_signal_cutoff = case_when(
        low_signal_rule == "blank_q99" ~ blank_q,
        low_signal_rule == "fixed" ~ low_signal_fixed,
        TRUE ~ NA_real_
      ),
      low_signal = !is.na(low_signal_cutoff) & mean_raw <= low_signal_cutoff,
      cv_primary = if_else(low_signal, cv_log, cv_raw),
      flag = case_when(
        n < 2 ~ "warn_missing_rep",
        is.na(cv_primary) ~ "warn_uncomputable",
        cv_primary >= cv_threshold ~ "fail_high_cv",
        TRUE ~ "pass"
      )
    ) %>%
    select(
      plate_id, run_date, plate_batch,
      sample, dilution, analyte,
      rep1, rep2,
      mean_raw, cv_raw,
      mean_log, cv_log,
      low_signal, low_signal_cutoff,
      cv_primary, flag
    )
  
  reps
}