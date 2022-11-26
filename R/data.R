#' Trust in scientific experts during the Covid-19 pandemic
#'
#' A dataset containing a measurement of "trust" in experts along with other
#' metrics collected through the Delphi Group at Carnegie Mellon University
#' U.S. COVID-19 Trends and Impact Survey, in partnership with Facebook. This
#' particular dataset is created from one of the public
#' [contingency tables](https://www.cmu.edu/delphi-web/surveys/monthly-rollup/),
#' specifically, the breakdown by state, age, gender, and race/ethnicity published
#' on 05 February 2022.
#'
#' @format A `data.frame` with 9759 rows and 8 columns
#' \describe{
#'   \item{`trust_experts`}{Real-valued. This is the average of
#'     `pct_trust_covid_info_*`
#'     where `*` is each of `doctors`, `experts`, `cdc`, and `govt_health`.}
#'   \item{`period`}{Factor. Start date of data collection period.
#'     There are 13 monthly periods}
#'   \item{`region`}{Factor. State abbreviation.}
#'   \item{`age`}{Factor. Self-reported age bucket.}
#'   \item{`gender`}{Factor. Self-reported gender.}
#'   \item{`raceethnicity`}{Factor. Self-reported race or ethnicity.}
#'   \item{`cli`}{Real-valued. This is the `wcli` indicator measuring the
#'     percent of circulating Covid-like illness in a particular region. See
#'     the [Delphi Epidata API](https://cmu-delphi.github.io/delphi-epidata/api/covidcast-signals/fb-survey.html#ili-and-cli-indicators)
#'     for a complete description.}
#'   \item{`hh_cmnty_cli`}{Real-valued. This is the `whh_cmnty_cli` indicator
#'     measuring the percent of people reporting illness in their local
#'     community and household.}
#' }
#' @source The U.S. COVID-19 Trends and Impact Survey.
#'
#'   The paper describing the survey:
#'
#'   Joshua A. Salomon, Alex Reinhart, Alyssa Bilinski, Eu Jing Chua,
#'   Wichada La Motte-Kerr, Minttu M. Rönn, Marissa Reitsma,
#'   Katherine Ann Morris, Sarah LaRocca, Tamar Farag, Frauke Kreuter,
#'   Roni Rosenfeld, and Ryan J. Tibshirani (2021). "The US COVID-19 Trends
#'   and Impact Survey: Continuous real-time measurement of COVID-19 symptoms,
#'   risks, protective behaviors, testing, and vaccination", Proceedings of the
#'   National Academy of Sciences 118 (51) e2111454118.
#'   \doi{10.1073/pnas.2111454118}.
#'
#'   [The Public Delphi US CTIS Documentation](https://cmu-delphi.github.io/delphi-epidata/symptom-survey/contingency-tables.html)
#'
#' @examples
#' \dontrun{
#' library(splines)
#' library(dplyr)
#' library(magrittr)
#' df <- 10
#'
#' trust_experts <- trust_experts %>%
#'   mutate(across(
#'     where(is.factor),
#'     ~ set_attr(.x, "contrasts", contr.sum(nlevels(.x), FALSE, TRUE))
#'   ))
#'
#' x <- Matrix::sparse.model.matrix(
#'     ~ 0 + region + age + gender + raceethnicity + period +
#'     bs(cli, df = df) + bs(hh_cmnty_cli, df = df),
#'     data = trust_experts, drop.unused.levels = TRUE)
#'
#' gr <- sapply(trust_experts, function(x) ifelse(is.factor(x), nlevels(x), NA))
#' gr <- rep(seq(ncol(trust_experts) - 1), times = c(gr[!is.na(gr)], df, df))
#' fit <- cv.sparsegl(x, trust_experts$trust_experts, gr)
#' }
"trust_experts"
