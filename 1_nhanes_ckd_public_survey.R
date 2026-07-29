# nhanes_ckd_public_survey.R
# Purpose:
#   Build a person-level analytic table from public NHANES files,
#   construct eGFR, urine ACR, and a binary CKD indicator,
#   and produce basic missingness and prevalence summaries.
#
# Notes:
#   This is a compact reproducible example. Full epidemiologic inference
#   would require survey weights, strata, and PSU design variables.

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(readr)
})

base_url <- "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018"
out_dir  <- "output"
dir.create(out_dir, showWarnings = FALSE)

read_nhanes_xpt <- function(file_stub) {
  url <- sprintf("%s/%s_J.XPT", base_url, file_stub)
  read_xpt(url)
}

# Public NHANES modules
demo   <- read_nhanes_xpt("DEMO")
biopro <- read_nhanes_xpt("BIOPRO")
albcr  <- read_nhanes_xpt("ALB_CR")
diq    <- read_nhanes_xpt("DIQ")
bpq    <- read_nhanes_xpt("BPQ")
bmx    <- read_nhanes_xpt("BMX")

analytic <- demo %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH1, WTMEC2YR, SDMVPSU, SDMVSTRA) %>%
  left_join(biopro %>% select(SEQN, LBXSCR), by = "SEQN") %>%
  left_join(albcr  %>% select(SEQN, URXUMA, URXUCR, URDACT), by = "SEQN") %>%
  left_join(diq    %>% select(SEQN, DIQ010), by = "SEQN") %>%
  left_join(bpq    %>% select(SEQN, BPQ020), by = "SEQN") %>%
  left_join(bmx    %>% select(SEQN, BMXBMI), by = "SEQN") %>%
  mutate(
    age = RIDAGEYR,
    sex = case_when(
      RIAGENDR == 1 ~ "Male",
      RIAGENDR == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    race_ethnicity = case_when(
      RIDRETH1 == 1 ~ "Mexican American",
      RIDRETH1 == 2 ~ "Other Hispanic",
      RIDRETH1 == 3 ~ "Non-Hispanic White",
      RIDRETH1 == 4 ~ "Non-Hispanic Black",
      RIDRETH1 == 5 ~ "Other / Multiracial",
      TRUE ~ NA_character_
    ),
    diabetes = case_when(
      DIQ010 == 1 ~ 1L,
      DIQ010 == 2 ~ 0L,
      TRUE ~ NA_integer_
    ),
    hypertension = case_when(
      BPQ020 == 1 ~ 1L,
      BPQ020 == 2 ~ 0L,
      TRUE ~ NA_integer_
    ),

    # 2021 CKD-EPI race-free creatinine equation
    female = if_else(sex == "Female", 1, 0, missing = 0),
    kappa  = if_else(female == 1, 0.7, 0.9),
    alpha  = if_else(female == 1, -0.241, -0.302),
    scr_k  = LBXSCR / kappa,

    egfr = 142 *
      pmin(scr_k, 1)^alpha *
      pmax(scr_k, 1)^(-1.200) *
      0.9938^age *
      1.012^female,

    # NHANES provides URDACT directly. If recomputing from raw units:
    # urine_acr = 100 * URXUMA / URXUCR
    urine_acr = URDACT,

    ckd_flag = case_when(
      !is.na(egfr) & !is.na(urine_acr) & (egfr < 60 | urine_acr >= 30) ~ 1L,
      !is.na(egfr) & !is.na(urine_acr) ~ 0L,
      TRUE ~ NA_integer_
    )
  )

missingness <- analytic %>%
  summarise(
    n = n(),
    serum_creatinine_missing = sum(is.na(LBXSCR)),
    urine_acr_missing        = sum(is.na(urine_acr)),
    egfr_missing             = sum(is.na(egfr)),
    ckd_flag_missing         = sum(is.na(ckd_flag))
  )

ckd_summary <- analytic %>%
  filter(!is.na(ckd_flag)) %>%
  summarise(
    n_complete = n(),
    crude_ckd_rate = mean(ckd_flag),
    diabetes_rate = mean(diabetes, na.rm = TRUE),
    hypertension_rate = mean(hypertension, na.rm = TRUE)
  )

write_csv(analytic,    file.path(out_dir, "nhanes_ckd_analytic.csv"))
write_csv(missingness, file.path(out_dir, "nhanes_ckd_missingness.csv"))
write_csv(ckd_summary, file.path(out_dir, "nhanes_ckd_summary.csv"))

print(missingness)
print(ckd_summary)
