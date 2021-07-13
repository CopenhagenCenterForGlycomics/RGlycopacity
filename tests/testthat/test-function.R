test_that("all main args must have non null values", {
  expect_error(main(data = NULL))
  expect_error(main(reference_data = NULL))
  expect_error(main(organism = NULL))
  expect_error(main(tissues = NULL))
  expect_error(main(quantile_cutoff = NULL))
  expect_error(main(reference_gene = NULL))
  expect_error(main(genes = NULL))
  expect_error(main(log_transform = NULL))
})

test_that("input matrix must be a matrix object", {
  expect_error(main(data = NA))
  expect_error(main(data = 'foo'))
  expect_error(main(data = 0))
  expect_error(main(data = c(0,1)))
  expect_error(main(data = data.frame(0)))
})
