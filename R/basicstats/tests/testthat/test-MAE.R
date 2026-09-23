test_that("MAE works", {
  obs <- c(1, 2, 3, 4, 5)
  pred <- c(1, 2, 3, 4, 5)
  expect_equal(MAE(obs, pred), 0)
})
