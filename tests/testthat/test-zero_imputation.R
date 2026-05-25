# Tests for zero_imputation

test_that("zero_imputation replaces zeros with ones", {
  RIBO <- matrix(
    c(0, 2, 0,
      4, 5, 6),
    nrow = 2, byrow = T,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )
  RNA <- matrix(
    c(7, 8, 9,
      0, 0.1, 3),
    nrow = 2, byrow = T,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )

  res <- zero_imputation(RIBO, RNA)

  expect_true(all(res$ribo >= 1))
  expect_true(all(res$rna  > 0))
  expect_equal(res$ribo[1, 1], 1)
  expect_equal(res$ribo[1, 3], 1)
  expect_equal(res$rna [2, 1], 0.1)
  expect_equal(res$rna [2, 2], 0.1)
})

test_that("zero_imputation preserves non-zero values and structure", {
  RIBO <- data.frame(
    s1 = c(0, 5),
    s2 = c(2, 0),
    row.names = c("g1", "g2")
  )
  RNA <- data.frame(
    s1 = c(3, 4),
    s2 = c(0, 6),
    row.names = c("g1", "g2")
  )

  res <- zero_imputation(RIBO, RNA)

  expect_identical(dim(res$ribo), dim(as.matrix(RIBO)))
  expect_identical(dim(res$rna),  dim(as.matrix(RNA)))
  expect_equal(res$ribo["g2", "s1"], 5)
  expect_equal(res$rna ["g2", "s2"], 6)
})

test_that("multiplicative imputation delegates to zCompositions", {
  testthat::skip_if_not_installed("zCompositions")

  RIBO <- matrix(
    c(0, 5, 1,
      2, 1, 0,
      3, 0, 4),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:3))
  )
  RNA <- matrix(
    c(1, 4, 0,
      0, 2, 5,
      3, 0, 1),
    nrow = 3, byrow = TRUE,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:3))
  )

  res <- zero_imputation(
    RIBO,
    RNA,
    method = "multiplicative",
    multiplicative_method = "GBM",
    output = "p-count"
  )

  expect_identical(dim(res$ribo), dim(RIBO))
  expect_identical(dim(res$rna),  dim(RNA))
  expect_equal(colnames(res$ribo), colnames(RIBO))
  expect_equal(rownames(res$ribo), rownames(RIBO))
  expect_true(all(res$ribo > 0))
  expect_true(all(res$rna  > 0))
})

test_that("simple zero_imputation caps replacement at 1 when minimum non-zero value is greater than 1", {
  RIBO <- matrix(
    c(0, 2, 5,
      3, 4, 6),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )
  RNA <- matrix(
    c(7, 0, 8,
      2, 3, 4),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )

  res <- zero_imputation(RIBO, RNA)

  # Smallest non-zero values are > 1, so zeros should still become 1
  expect_equal(res$ribo["g1", "s1"], 1)
  expect_equal(res$rna["g1", "s2"], 1)

  # Non-zero values should stay unchanged
  expect_equal(res$ribo["g1", "s2"], 2)
  expect_equal(res$ribo["g2", "s3"], 6)
  expect_equal(res$rna["g2", "s1"], 2)
  expect_equal(res$rna["g1", "s3"], 8)

  # Everything should be positive after imputation
  expect_true(all(res$ribo > 0))
  expect_true(all(res$rna > 0))
})
#
