# run_render_folder.R
library(quarto)

folder <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(folder)) folder <- "data/raw"

files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)

dir.create("outputs/reports", recursive = TRUE, showWarnings = FALSE)

for (f in files) {
  out_html <- file.path("outputs/reports", paste0(gsub("\\.csv$", "", basename(f)), "_QC.html"))
  quarto_render(
    input = "QCReport_HEV.qmd",
    execute_params = list(
      plate_file = f,
      value_datatype = "Net MFI",
      cv_threshold = 0.30,
      blank_quantile = 0.99,
      min_total_events = 50
    ),
    output_file = out_html
  )
  message("Rendered: ", out_html)
}