# nhanes_weighted_correlation_table.R
# Purpose:
#   Create an ordinal CKD severity variable and compute unweighted
#   and sample-weighted correlations with candidate variables.
#
#   Weighted correlations here use MEC weights only. Full complex-survey
#   standard errors would additionally use strata and PSU variables.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

infile  <- "output/nhanes_ckd_analytic.csv"
out_dir <- "output"

dat <- read_csv(infile, show_col_types = FALSE)

dat <- dat %>%
  mutate(
    ckd_stage = case_when(
      is.na(egfr) | is.na(urine_acr) ~ NA_real_,
      egfr >= 60 & urine_acr < 30   ~ 0,
      egfr >= 90 & urine_acr >= 30  ~ 1,
      egfr >= 60 & egfr < 90 & urine_acr >= 30 ~ 2,
      egfr >= 45 & egfr < 60        ~ 3,
      egfr >= 30 & egfr < 45        ~ 4,
      egfr >= 15 & egfr < 30        ~ 5,
      egfr < 15                     ~ 6
    ),
    male = if_else(sex == "Male", 1, 0, missing = NA_real_)
  )

weighted_cor <- function(x, y, w) {
  keep <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- x[keep]
  y <- y[keep]
  w <- w[keep]

  if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
    return(NA_real_)
  }

  wx <- sum(w * x) / sum(w)
  wy <- sum(w * y) / sum(w)

  cov_xy <- sum(w * (x - wx) * (y - wy)) / sum(w)
  var_x  <- sum(w * (x - wx)^2) / sum(w)
  var_y  <- sum(w * (y - wy)^2) / sum(w)

  cov_xy / sqrt(var_x * var_y)
}

candidate_vars <- c(
  "age", "male", "BMXBMI", "diabetes", "hypertension",
  "LBXSCR", "urine_acr", "egfr"
)

candidate_vars <- candidate_vars[candidate_vars %in% names(dat)]

cor_table <- lapply(candidate_vars, function(v) {
  x <- dat[[v]]
  y <- dat$ckd_stage

  complete <- is.finite(x) & is.finite(y)

  tibble(
    variable = v,
    n_complete = sum(complete),
    pearson_r = if (sum(complete) >= 3) cor(x[complete], y[complete]) else NA_real_,
    weighted_r = weighted_cor(x, y, dat$WTMEC2YR)
  )
}) %>%
  bind_rows() %>%
  mutate(abs_pearson_r = abs(pearson_r)) %>%
  arrange(desc(abs_pearson_r))

write_csv(cor_table, file.path(out_dir, "nhanes_ckd_correlation_table.csv"))

print(cor_table)
