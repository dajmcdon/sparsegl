library(tidyverse)
## Downloaded on 10 February 2022
covidcast_url <- "https://www.cmu.edu/delphi-web/surveys/monthly-rollup/monthly_state_age_gender_raceethnicity.csv.gz"
symp <- read_csv(covidcast_url)
symp <- symp %>%
  select(period_start, region, age, gender, raceethnicity, starts_with("val_"))
sm <- symp %>%
  select(period_start, region, age, gender, raceethnicity,
         contains("val_pct_trust_covid_info"),
         val_pct_cli,
         val_pct_hh_cmnty_cli,
         val_pct_wearing_mask_5d, val_pct_wearing_mask_7d) %>%
  mutate(period = lubridate::ymd(period_start)) %>%
  filter(period > lubridate::ymd("2021-05-01")) %>% # remove pre-survey period
  select(-period_start) %>%
  mutate(
    age = str_c("age_", replace_na(age, "NotReported")),
    gender = str_c("gender_", replace_na(gender, "NotReported")),
    raceethnicity = str_c("race_", replace_na(raceethnicity, "NotReported"))
  ) %>%
  rename_with(~str_remove(.x, "val_pct_"), starts_with("val_pct"))

# Process to x, y ---------------------------------------------------------

one_hot <- function(.data, ...) {
  ex <- rlang::expr(c(...))
  pos <- names(tidyselect::eval_select(ex, .data)) # gives names of columns
  .data <- mutate(.data, .id = 1:n()) # add an id var just to be sure
  for (p in pos) {
    var <- sym(p)
    .data <- mutate(.data, temp = 1) %>%
      pivot_wider(names_from = !!var, values_from = temp, values_fill = 0)
  }
  return(select(.data, -.id))
}

my_splines <- function(x, .name, df = 10) {
  as_tibble(splines::bs(x, df=df)) %>%
    setNames(str_c(.name, 1:df))
}


sm <- sm %>%
  one_hot(period, region, age, gender, raceethnicity) %>%
  select(!(trust_covid_info_politicians:trust_covid_info_religious)) %>%
  mutate(my_splines(cli, "cli_"),
         my_splines(hh_cmnty_cli, "cmnty_cli_"),
         .keep = "unused") %>%
  rowwise() %>%
  mutate(
    trust_experts = mean(c_across(starts_with("trust_covid")), na.rm=TRUE),
    masking = mean(c_across(starts_with("wearing_mask")), na.rm=TRUE)) %>%
  select(! starts_with("trust_covid") & ! starts_with("wearing_mask")) %>%
  select(! masking)



cc <- complete.cases(sm_onehot)
y <- sm_onehot$trust_experts[cc]
x <- Matrix::Matrix(
  as.matrix(sm_onehot %>% select(!(trust_experts)))[cc, ],
  sparse = TRUE)

trust_experts <- cbind(y, x)

usethis::use_data(trust_experts, overwrite = TRUE)
