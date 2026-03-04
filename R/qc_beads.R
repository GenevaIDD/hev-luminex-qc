# R/qc_beads.R
library(dplyr)

qc_low_beads <- function(df, bead_threshold = 30) {
  df %>%
    filter(!is.na(bead_count)) %>%
    mutate(low_beads = bead_count < bead_threshold) %>%
    summarise(
      n_points = n(),
      n_low = sum(low_beads, na.rm = TRUE),
      pct_low = mean(low_beads, na.rm = TRUE),
      .by = c(plate_id)
    ) %>%
    left_join(
      df %>%
        filter(!is.na(bead_count), bead_count < bead_threshold) %>%
        distinct(plate_id, well, sample, dilution, sample_type, analyte, bead_count) %>%
        arrange(bead_count),
      by = "plate_id"
    )
}