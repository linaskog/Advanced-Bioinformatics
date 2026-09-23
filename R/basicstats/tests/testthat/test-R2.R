test_that("R2 works", {
  obs <- c(1, 2, 3, 4, 5)
  pred <- c(1, 2, 3, 4, 5)
  except_equal(R2(obs, pred), 1)
})
