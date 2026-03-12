# R/standard_4pl.R
library(dplyr)

# 4PL function: y = d + (a-d)/(1 + (x/c)^b)
f_4pl <- function(x, a, b, c, d) {
  d + (a - d) / (1 + (x / c)^b)
}

# Fit one analyte's 4PL on one plate
fit_4pl_one <- function(dat) {
  dat <- dat %>%
    filter(!is.na(denom), denom > 0, !is.na(value)) %>%
    arrange(denom)

  if (dplyr::n_distinct(dat$denom) < 4) return(NULL)

  # Starting values for decreasing dose-response curve:
  #   a = upper asymptote (high conc / low denom → high MFI)
  #   d = lower asymptote (low conc / high denom → low MFI, near 0 after blank subtraction)
  #   c = inflection point (midpoint on denom scale)
  #   b = slope (positive for decreasing curve with this parameterisation)
  a0 <- max(dat$value, na.rm = TRUE)
  d0 <- min(dat$value, na.rm = TRUE)
  c0 <- stats::median(dat$denom, na.rm = TRUE)
  b0 <- 1

  # Weights proportional to 1/y (not 1/y²).  The 1/y² scheme
  # over-aggressively downweights the highest-concentration standards,
  # causing the fit to miss the upper asymptote for analytes like
  # ORF2 Bohm and ORF2-E2 Bohm.  Using 1/y corresponds to a Poisson-like
  # variance model (var ∝ mean), which is appropriate for photon-counting-
  # based Luminex measurements and still upweights low-signal points.
  w <- 1 / pmax(dat$value, 1)

  # Constrain with port algorithm
  tryCatch(
    stats::nls(
      value ~ f_4pl(denom, a, b, c, d),
      data = dat,
      start = list(a = a0, b = b0, c = c0, d = d0),
      weights = w,
      algorithm = "port",
      lower = c(a = 0, b = 0.01, c = min(dat$denom) * 0.1, d = 0),
      upper = c(a = Inf, b = 10,   c = max(dat$denom) * 10,  d = Inf),
      control = nls.control(maxiter = 200, warnOnly = TRUE)
    ),
    error = function(e) NULL
  )
}

# Invert 4PL: given y, return x (denom-equivalent "relative concentration")
invert_4pl <- function(y, a, b, c, d) {
  eps <- 1e-6
  
  lo <- pmin(a, d) + eps
  hi <- pmax(a, d) - eps
  y2 <- pmin(pmax(y, lo), hi)
  
  inner <- (a - d) / (y2 - d) - 1
  
  # inner should be >= 0 for valid inversion
  inner[inner < 0] <- NA_real_
  
  c * (inner^(1 / b))
}


fit_4pl_all <- function(std_df) {
  # expects: plate_id, analyte, denom, value, is_blank
  std_only <- std_df %>%
    filter(!is_blank, !is.na(denom))

  # Observed MFI range of the standards (for LLOQ/ULOQ flagging)
  std_range <- std_only %>%
    group_by(plate_id, analyte) %>%
    summarise(
      std_mfi_min = min(value, na.rm = TRUE),
      std_mfi_max = max(value, na.rm = TRUE),
      .groups = "drop"
    )

  fits <- std_only %>%
    group_by(plate_id, analyte) %>%
    group_modify(~{
      mod <- fit_4pl_one(.x)
      if (is.null(mod)) {
        return(tibble(
          a = NA_real_, b = NA_real_, c = NA_real_, d = NA_real_,
          fit_ok = FALSE,
          fit_reason = "nls_failed"
        ))
      }
      co <- coef(mod)
      tibble(
        a = unname(co["a"]),
        b = unname(co["b"]),
        c = unname(co["c"]),
        d = unname(co["d"]),
        fit_ok = TRUE,
        fit_reason = NA_character_
      )
    }) %>%
    ungroup() %>%
    left_join(std_range, by = c("plate_id", "analyte"))

  fits
}

normalize_df_to_4pl <- function(df, fits) {
  # df expects: plate_id, analyte, value (Net MFI)
  # fits expects: plate_id, analyte, a, b, c, d, fit_ok, std_mfi_min, std_mfi_max
  #
  # Adds columns:
  #   norm_denom  – inverted 4PL concentration (NA when out of range or fit failed)
  #   quant_flag  – "ok", ">ULOQ", "<LLOQ", or NA (fit failed / missing value)
  df %>%
    left_join(fits, by = c("plate_id", "analyte")) %>%
    mutate(
      quant_flag = case_when(
        !fit_ok | is.na(value)          ~ NA_character_,
        value > std_mfi_max             ~ ">ULOQ",
        value < std_mfi_min             ~ "<LLOQ",
        TRUE                            ~ "ok"
      ),
      norm_denom = if_else(
        fit_ok & !is.na(value) & quant_flag == "ok",
        invert_4pl(value, a, b, c, d),
        NA_real_
      )
    )
}