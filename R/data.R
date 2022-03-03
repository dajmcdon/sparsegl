#' Trust in scientific experts during the Covid-19 pandemic
#'
#' A dataset containing a measurement of "trust" in experts along with other
#' metrics collected through the Delphi Group at Carnegie Mellon University
#' U.S. COVID-19 Trends and Impact Survey, in partnership with Facebook. This
#' particular dataset is created from one of the public [contingency tables](https://www.cmu.edu/delphi-web/surveys/monthly-rollup/),
#' specifically, the breakdown by state, age, gender, and race/ethnicity published
#' on 05 February 2022.
#'
#' @format A sparse [Matrix::sparseMatrix()] with 3775 rows, 96 columns, and
#'   51738 non-zero entries
#' \describe{
#'   \item{`y`}{Real-valued response. This is the average of `pct_trust_covid_info_*`
#'     where `*` is each of `doctors`, `experts`, `cdc`, and `govt_health`.}
#'   \item{`yyyy-mm-01`}{`0-1`-valued predictor. Start date of data collection period.
#'     There are 8 monthly periods}
#'   \item{`AK`-`WY`}{`0-1`-valued predictor. State abbreviation.}
#'   \item{`age_*`}{`0-1`-valued predictor. Self-reported age bucket.}
#'   \item{`gender_*`}{`0-1`-valued predictor. Self-reported gender.}
#'   \item{`race_*`}{`0-1`-valued predictor. Self-reported race.}
#'   \item{`cli_*`}{Real-valued predictor. `pct_cli` expanded in a B-spline
#'     basis with 10 degrees of freedom.}
#'   \item{`cmnty_cli_*`}{Real-valued predictor. `pct_hh_cmnty_cli` expanded
#'     in a B-spline basis with 10 degrees of freedom.}
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
"trust_experts"
