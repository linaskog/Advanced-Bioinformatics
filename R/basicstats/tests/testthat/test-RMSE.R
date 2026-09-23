test_that("RMSE works", {
  obs <- c(1, 2, 3, 4, 5)
  pred <- c(1, 2, 3, 4, 5)
  expect_equal(RMSE(obs, pred), 0)
})
