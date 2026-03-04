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
  
  # Robust-ish starting values
  a0 <- min(dat$value, na.rm = TRUE)
  d0 <- max(dat$value, na.rm = TRUE)
  c0 <- stats::median(dat$denom, na.rm = TRUE)
  b0 <- 1
  
  # Constrain with port algorithm
  tryCatch(
    stats::nls(
      value ~ f_4pl(denom, a, b, c, d),
      data = dat,
      start = list(a = a0, b = b0, c = c0, d = d0),
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
  
  fits <- std_only %>%
    group_by(plate_id, analyte) %>%
    group_modify(~{
      mod <- fit_4pl_one(.x)
      if (is.null(mod)) {
        return(tibble(
          a = NA_real_, b = NA_real_, c = NA_real_, d = NA_real_,
          fit_ok = FALSE,
          fit_reason = "nls_failed"
        ))      }
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
    ungroup()
  
  fits
}

normalize_df_to_4pl <- function(df, fits) {
  # df expects: plate_id, analyte, value (Net MFI)
  df %>%
    left_join(fits, by = c("plate_id","analyte")) %>%
    mutate(
      norm_denom = if_else(
        fit_ok & !is.na(value),
        invert_4pl(value, a, b, c, d),
        NA_real_
      )
    )
}