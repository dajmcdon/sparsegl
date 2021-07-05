test_that("group norms are correct", {
  asparse <- .05
  gr <- c(1,1,1,1,2,2,3,3,3,3,3)
  x <- c(-5):5
  expect_equal(grouped_one_norm(x, gr), c(14, 1, 15))
  expect_equal(grouped_two_norm(x, gr),
               c(two_norm(x[1:4]), two_norm(x[5:6]), two_norm(x[7:11])))
  expect_equal(sp_group_norm(x, gr, asparse),
               asparse*one_norm(x) + (1 - asparse)*gr_two_norm(x, gr))
  asparse <- .3
  expect_equal(sp_group_norm(x, gr, asparse),
               asparse*one_norm(x) + (1 - asparse)*gr_two_norm(x, gr))
  expect_error(sp_group_norm(x, 1, 0.05))
  expect_error(sp_group_norm(x, gr, 3))
  expect_error(sp_group_norm(x, gr, -1))
})
