# R/qc_replicates.R
library(dplyr)
library(tidyr)

compute_replicate_variability <- function(df,
                                          cv_threshold = 0.30,
                                          low_mfi_threshold = 50) {

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
      .groups = "drop"
    ) %>%
    mutate(
      low_signal = mean_raw < low_mfi_threshold,
      flag = case_when(
        n < 2                  ~ "warn_missing_rep",
        is.na(cv_raw)          ~ "warn_uncomputable",
        low_signal             ~ "low_signal",
        cv_raw >= cv_threshold ~ "fail_high_cv",
        TRUE                   ~ "pass"
      )
    ) %>%
    select(
      plate_id, run_date, plate_batch,
      sample, dilution, analyte,
      rep1, rep2,
      mean_raw, cv_raw,
      low_signal, flag
    )

  reps
}