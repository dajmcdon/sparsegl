test_that("family check returns correct results", {


  expect_identical(check_family("g"), list(check = "char", family = "g"))
  expect_identical(check_family(gaussian), list(check = "fam", family = gaussian()))
  fake_family <- gaussian()
  bad_class <- fake_family
  class(bad_class) <- "list"
  expect_identical(check_family(bad_class), list(check = "warn", family = bad_class))
  bad_class$linkinv <- NULL
  expect_identical(check_family(bad_class), list(check = "err", family = bad_class))
  expect_identical(
    check_family(binomial(link = "probit")),
    list(check = "fam", family = binomial(link = "probit"))
  )
  expect_error(check_family(mean))
  f <- function() 1:5
  expect_identical(check_family(f()), list(check = "err", family = f()))
  f <- function() list(a = "a")
  expect_identical(check_family(f()), list(check = "err", family = f()))
})
