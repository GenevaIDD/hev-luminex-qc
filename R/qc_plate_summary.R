# R/qc_plate_summary.R
library(dplyr)

plate_summary <- function(df,
                          min_total_events = 50,
                          blank_quantile = 0.99) {
  
  # Events
  events_tbl <- df %>%
    filter(!is.na(total_events)) %>%
    distinct(plate_id, well, total_events) %>%
    summarise(
      min_events = min(total_events, na.rm = TRUE),
      pct_low_events = mean(total_events < min_total_events, na.rm = TRUE),
      .by = plate_id
    )
  
  # Blank quantiles per analyte
  blank_tbl <- df %>%
    filter(sample_type == "blank") %>%
    group_by(plate_id, analyte) %>%
    summarise(
      blank_q99 = quantile(value, blank_quantile, na.rm = TRUE),
      blank_median = median(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Specimen counts
  counts <- 
    df %>%
    distinct(plate_id, well_type, sample, dilution) %>%
    group_by(well_type,dilution) %>%
    summarise(n = n())
  
  analytes <- 
    df %>%
    distinct(plate_id, Location, analyte) %>%
    group_by(analyte) %>%
    summarise(n = n())
  
  
  list(events = events_tbl, blanks = blank_tbl, counts = counts, analytes=analytes)
}