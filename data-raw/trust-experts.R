library(tidyverse)
## Downloaded on 18 November 2022
covidcast_url <- "https://www.cmu.edu/delphi-web/surveys/monthly-rollup/monthly_state_all_indicators_age_gender_raceethnicity.csv.gz"
symp <- read_csv(covidcast_url)
symp <- symp |>
  select(period_start, region, age, gender, raceethnicity, starts_with("val_"))
trust_experts <- symp |>
  select(
    period_start, region, age, gender, raceethnicity,
    contains("val_pct_trust_covid_info"),
    val_pct_cli,
    val_pct_hh_cmnty_cli,
    val_pct_wearing_mask_5d, val_pct_wearing_mask_7d
  ) |>
  mutate(period = lubridate::ymd(period_start)) |>
  filter(period > lubridate::ymd("2021-05-01")) |> # remove pre-survey period
  select(-period_start) |>
  mutate(
    age = str_c(replace_na(age, "NotReported")),
    gender = str_c(replace_na(gender, "NotReported")),
    raceethnicity = str_c(replace_na(raceethnicity, "NotReported"))
  ) |>
  rename_with(~ str_remove(.x, "val_pct_"), starts_with("val_pct"))



trust_experts <- trust_experts |>
  select(!(trust_covid_info_politicians:trust_covid_info_religious)) |>
  mutate(
    region = as.factor(region),
    age = as.factor(age),
    gender = as.factor(gender),
    raceethnicity = as.factor(raceethnicity),
    period = as.factor(period)
  ) |>
  rowwise() |>
  mutate(
    trust_experts = mean(c_across(starts_with("trust_covid")), na.rm = TRUE),
    masking = mean(c_across(starts_with("wearing_mask")), na.rm = TRUE)
  ) |>
  select(!starts_with("trust_covid") & !starts_with("wearing_mask")) |>
  select(!masking) |>
  ungroup()

cc <- complete.cases(trust_experts)
trust_experts <- as.data.frame(trust_experts[cc, ])

usethis::use_data(trust_experts, overwrite = TRUE)
