# R/subtract_blanks.R
#
# Blank subtraction for true Net MFI.
# The xPONENT "Net MFI" block is actually just Median MFI because blanks
# were not identified in the instrument. This function subtracts the mean
# of the blank wells (per plate/analyte) from all values to produce a
# true background-corrected Net MFI. Negative values are floored at 0.

library(dplyr)

subtract_plate_blanks <- function(df) {
  blank_means <- df %>%
    filter(sample_type == "blank") %>%
    group_by(plate_id, analyte) %>%
    summarise(blank_mean = mean(value, na.rm = TRUE), .groups = "drop")

  df %>%
    left_join(blank_means, by = c("plate_id", "analyte")) %>%
    mutate(
      value = pmax(value - blank_mean, 0)
    ) %>%
    select(-blank_mean)
}
