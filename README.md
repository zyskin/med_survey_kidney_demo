# med_survey_kidney_demo

These code samples use public NHANES survey data as a compact example of a reproducible survey-data workflow: merge respondent-level files, construct an analytic table, retain survey-design variables, check missingness, compute derived indicators, and run basic statistical/ML analysis.

NHANES is a complex U.S. health survey with interviews, physical exams, and lab measurements. In this example, the target variable is a simple CKD indicator derived from serum creatinine/eGFR and urine albumin-creatinine ratio. The health topic is only a demonstration; the transferable part is the survey-data processing and modelling workflow.

**Code 1 — R: data ingestion and analytic-file construction**
Downloads/reads NHANES modules, merges them by respondent ID, recodes variables, computes eGFR/ACR/CKD flag, and writes a cleaned person-level file. It retains `WTMEC2YR`, `SDMVPSU`, and `SDMVSTRA`, but the printed summaries are crude/unweighted.

**Code 2 — R: weighted correlation table**
Uses the cleaned file from Code 1 and computes correlations between an ordinal CKD severity variable and candidate predictors. This is the code that actually uses weights: `WTMEC2YR` is passed into the weighted-correlation function. It uses weighted means, covariance, and variances. It does not implement full complex-survey variance estimation with PSU/strata.

**Code 3 — Python: ML classification pipeline**
Reads the cleaned file, defines a binary CKD outcome, performs train/test splitting, imputation, scaling, logistic regression, random forest fitting, and test-set evaluation. This sample is unweighted; it demonstrates reproducible modelling workflow rather than formal population estimation.

`WTMEC2YR` is the NHANES two-year MEC examination weight, taken from the NHANES demographic file. Since the example uses exam/lab variables, it is the relevant weight to retain and use for weighted estimates. For formal NHANES inference, better to use `WTMEC2YR` together with `SDMVPSU` and `SDMVSTRA`.
