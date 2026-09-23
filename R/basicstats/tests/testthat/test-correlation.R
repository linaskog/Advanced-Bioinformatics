test_that("correlation works", {
  obs <- c(1, 2, 3)
  pred <- c(2, 4, 6)
  expect_equal(correlation(obs, pred), 1)
})
