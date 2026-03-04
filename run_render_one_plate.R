# run_render_one_plate.R
# Render a single HEV Luminex plate QC report (Quarto)
# Workaround: some Quarto builds don't allow --output to include a path.

library(quarto)
library(tools)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_render_one_plate.R <path_to_plate_csv>")
}

plate_file <- args[1]
if (!file.exists(plate_file)) stop("Plate file does not exist: ", plate_file)

dir.create("outputs/reports", recursive = TRUE, showWarnings = FALSE)

# Desired final report name
final_html <- file.path(
  "outputs/reports",
  paste0(file_path_sans_ext(basename(plate_file)), "_QC.html")
)

# Render with output filename ONLY (no path)
tmp_name <- paste0(file_path_sans_ext(basename(plate_file)), "_QC.html")

# Render in current working directory
quarto_render(
  input = "QCReport_HEV.qmd",
  execute_params = list(
    plate_file = plate_file,
    value_datatype = "Net MFI",
    cv_threshold = 0.30,
    blank_quantile = 0.99,
    min_total_events = 50
  ),
  output_file = tmp_name
)

# Move to outputs/reports
if (!file.exists(tmp_name)) stop("Expected output not found: ", tmp_name)
file.rename(tmp_name, final_html)

cat("QC report written to:", final_html, "\n")