# R/read_xponent_plate.R
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)

source("R/utils.R")

.extract_header_field <- function(lines, key) {
  hit <- lines[str_detect(lines, paste0("^", fixed(key), ","))]
  if (length(hit) == 0) return(NA_character_)
  str_split(hit[1], ",", simplify = TRUE)[1,2] |> str_trim()
}

.find_datatype_start <- function(lines, datatype) {
  which(str_detect(lines, paste0("^DataType:,", fixed(datatype), "\\b")))
}

.read_block_as_df <- function(lines, start_line) {
  # start_line points to "DataType:,..."
  header <- readr::read_csv(I(lines[start_line + 1]),
                            col_names = FALSE, show_col_types = FALSE) |>
    unlist(use.names = FALSE) |>
    as.character()
  
  i <- start_line + 2
  out <- list()
  
  while (i <= length(lines)) {
    if (lines[i] == "" || stringr::str_detect(lines[i], "^DataType:,")) break
    if (stringr::str_detect(lines[i], "^,+$")) { i <- i + 1; next }
    
    out[[length(out) + 1]] <- readr::read_csv(
      I(lines[i]),
      col_names = FALSE,
      show_col_types = FALSE
    ) |>
      unlist(use.names = FALSE)
    
    i <- i + 1
  }
  
  if (length(out) == 0) return(tibble::tibble())
  
  # Build rows with placeholder names V1..Vn (so tibble is happy)
  df <- dplyr::bind_rows(lapply(out, function(x) {
    tmp <- as.list(x)
    names(tmp) <- paste0("V", seq_along(tmp))
    tibble::as_tibble_row(tmp)
  }))
  
  # Ensure header length matches df columns
  if (length(header) < ncol(df)) header <- c(header, rep("", ncol(df) - length(header)))
  header <- header[seq_len(ncol(df))]
  
  # Repair blank/duplicate header names
  # - Make blanks into X1, X2, ...
  # - Make duplicates unique
  header <- make.unique(ifelse(is.na(header) | header == "", "X", header), sep = "_")
  
  names(df) <- header
  df
}
.parse_well <- function(location) {
  # "1(1,A1)" -> "A1"
  str_match(location, "\\(([A-H]\\d{1,2})\\)$")[,2]
}

.classify_sample <- function(sample) {
  case_when(
    str_detect(sample, regex("^blank$", ignore_case = TRUE)) ~ "blank",
    str_detect(sample, regex("^std\\b", ignore_case = TRUE)) ~ "standard",
    TRUE ~ "specimen"
  )
}

.extract_dilution <- function(sample, dilution_col) {
  d1 <- na_if(str_trim(as.character(dilution_col)), "")
  d2 <- str_match(sample, "(\\d+\\s*/\\s*\\d+)")[,2] |> str_replace_all("\\s+", "")
  coalesce(d1, d2)
}

classify_well_type <- function(df, sample_col = "sample_name") {
  
  s <- str_to_lower(as.character(df[[sample_col]]))
  
  df %>%
    mutate(
      well_type = case_when(
        str_detect(s, "^std\\b|\\bstandard\\b")        ~ "standard",
        str_detect(s, "\\bpos\\b|positive control")   ~ "control_pos",
        str_detect(s, "\\bneg\\b|negative control")   ~ "control_neg",
        str_detect(s, "\\bblank\\b|background|bkg")   ~ "blank",
        TRUE                                           ~ "specimen"
      )
    )
}

read_xponent_plate <- function(path,
                               value_datatype = "Net MFI",
                               id_datatype = "Median") {
  
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  
  plate_id <- guess_plate_id(path)
  plate_batch <- .extract_header_field(lines, "Batch")
  run_date <- .extract_header_field(lines, "Date")
  panel <- .extract_header_field(lines, "PanelName")
  
  # 1) ID mapping from Median (has Sample IDs + dilution col in your files)
  s_id <- .find_datatype_start(lines, id_datatype)
  stopifnot(length(s_id) >= 1)
  id_df <- .read_block_as_df(lines, s_id[1])
  
  # In your plates the 3rd column header is blank/placeholder; grab it by position safely
  dil_col <- names(id_df)[3]
  
  id_map <- id_df |>
    transmute(
      location_raw = .data$Location,
      well = .parse_well(.data$Location),
      sample = as.character(.data$Sample),
      dilution_raw = .data[[dil_col]],
      dilution = .extract_dilution(sample, dilution_raw),
      sample_type = .classify_sample(sample)
    )
  
  # 2) Values from chosen DataType, then join back on Location
  s_val <- .find_datatype_start(lines, value_datatype)
  stopifnot(length(s_val) >= 1)
  val_df <- .read_block_as_df(lines, s_val[1])
  
  # Identify Total Events column if present
  has_events <- "Total Events" %in% names(val_df)
  
  val_df2 <- val_df |>
    mutate(
      well = .parse_well(.data$Location),
      total_events = if (has_events) safe_numeric(.data$`Total Events`) else NA_real_
    )
  
  # Long format: analyte columns are everything except core columns
  core_cols <- c("Location", "Sample", "Total Events", "well", "total_events")
  keep_core <- intersect(core_cols, names(val_df2))
  
  val_long <- val_df2 |>
    pivot_longer(
      cols = setdiff(names(val_df2), keep_core),
      names_to = "analyte",
      values_to = "value"
    ) |>
    mutate(value = safe_numeric(value)) |>
    #  DROP NON-ANALYTE PLACEHOLDERS
    filter(
      !stringr::str_detect(analyte, "^X(_\\d+)?$"),
      analyte != "",
      !is.na(analyte)
    )  

  
  # --- ALSO read bead counts (DataType: Count) ---
  s_cnt <- .find_datatype_start(lines, "Count")
  
  bead_long <- tibble::tibble(Location = character(), well = character(),
                              analyte = character(), bead_count = numeric())
  
  ## getting rid of placeholder analytes like "X", "X_1", etc.
  bead_long <- bead_long |>
    dplyr::filter(!stringr::str_detect(analyte, "^X(_\\d+)?$"))
  
  if (length(s_cnt) >= 1) {
    cnt_df <- .read_block_as_df(lines, s_cnt[1])
    
    has_events_cnt <- "Total Events" %in% names(cnt_df)
    
    cnt_df2 <- cnt_df |>
      mutate(
        well = .parse_well(.data$Location)
      )
    
    core_cnt <- c("Location", "Sample", "Total Events", "well")
    keep_core_cnt <- intersect(core_cnt, names(cnt_df2))
    
    bead_long <- cnt_df2 |>
      pivot_longer(
        cols = setdiff(names(cnt_df2), keep_core_cnt),
        names_to = "analyte",
        values_to = "bead_count"
      ) |>
      mutate(bead_count = safe_numeric(bead_count)) |>
      select(Location, well, analyte, bead_count)
  }
  
  # join bead counts onto Net MFI long data
  val_long <- val_long |>
    left_join(bead_long, by = c("Location", "well", "analyte"))
  
  out <- val_long |>
    left_join(id_map, by = c("Location" = "location_raw", "well" = "well")) |>
    mutate(
      plate_file = basename(path),
      plate_id = plate_id,
      plate_batch = plate_batch,
      run_date = run_date,
      panel = panel,
      value_datatype = value_datatype
    ) |>
    classify_well_type(sample_col = "sample") |>
    group_by(plate_id, sample_type, sample, dilution, analyte) |>
    mutate(rep_index = row_number()) |>
    ungroup()
  
  out
}