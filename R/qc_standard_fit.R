# R/qc_standard_fit.R

library(dplyr)
library(stringr)
library(purrr)
library(nls.multstart)

# ---- Fit 4PL per plate × analyte ----
fit_4pl <- function(std_df) {
  
  # Expect columns: analyte, denom, value
  std_df <- std_df %>% filter(!is.na(denom), value > 0)
  
  if (n_distinct(std_df$denom) < 4) {
    return(NULL)  # not enough points for 4PL
  }
  
  tryCatch(
    nls_multstart(
      value ~ d + (a - d) / (1 + (denom / c)^b),
      data = std_df,
      iter = 500,
      start_lower = c(a = 0, d = max(std_df$value) * 0.5, c = min(std_df$denom), b = 0.5),
      start_upper = c(a = min(std_df$value), d = max(std_df$value) * 2, c = max(std_df$denom), b = 5),
      supp_errors = "Y",
      na.action = na.omit,
      control = nls.control(maxiter = 200)
    ),
    error = function(e) NULL
  )
}

# ---- Invert 4PL: MFI -> relative concentration ----
invert_4pl <- function(y, coef) {
  a <- coef["a"]
  b <- coef["b"]
  c <- coef["c"]
  d <- coef["d"]
  
  # Guardrails
  y <- pmin(pmax(y, a + 1e-6), d - 1e-6)
  
  c * (( (a - d) / (y - d) - 1 )^(1 / b))
}

# ---- Apply per plate × analyte ----
normalize_to_standard <- function(df, std_df) {
  
  fits <- std_df %>%
    group_by(plate_id, analyte) %>%
    group_map(~ {
      mod <- fit_4pl(.x)
      if (is.null(mod)) return(NULL)
      
      tibble(
        plate_id = .y$plate_id,
        analyte = .y$analyte,
        coef = list(coef(mod))
      )
    }) %>%
    bind_rows()
  
  if (nrow(fits) == 0) return(df)
  
  df %>%
    left_join(fits, by = c("plate_id", "analyte")) %>%
    mutate(
      norm_conc = if_else(
        !is.na(coef),
        map2_dbl(value, coef, invert_4pl),
        NA_real_
      )
    ) %>%
    select(-coef)
}