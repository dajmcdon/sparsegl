
check_family <- function(family) {
  if (is.character(family)) return(list(check = "char", family = family))
  if (is.function(family)) family <- family()
  has_family_class <- inherits(family, "family")
  if (!is.list(family)) return(list(check = "err", family = family))
  has_vfun <- is.function(family$variance)
  has_invlink <- is.function(family$linkinv)

  if (!(has_vfun && has_invlink)) check <- "err"
  else if (!has_family_class) check <- "warn"
  else check <- "fam"
  list(check = check, family = family)
}

validate_family <- function(family) {
  check <- check_family(family)
  if (check$check == "warn") {
    cli::cli_warn(c(
      "`family` does not have class {.cls family}, but appears to contain",
      i = "the required functions {.field variance} and {.field linkinv}.",
      i = "Attempting to estimate sparse group lasso with IRLS."
    ))
  }
  if (check$check == "err") {
    cli::cli_abort("`family` is not to be a valid family object. See `?family`.")
  }
  invisible(check)
}
