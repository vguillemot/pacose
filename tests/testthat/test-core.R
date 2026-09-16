test_that("partial correlations convert to precision and back", {
  partial <- matrix(c(1, 0.2, -0.3,
                      0.2, 1, 0.4,
                      -0.3, 0.4, 1), nrow = 3, byrow = TRUE)
  variances <- c(2, 3, 4)

  precision <- pcor2invcov(partial, variances)

  expect_equal(invcov2pcor(precision), partial, tolerance = 1e-12)
})

test_that("oppdiag changes only the diagonal sign", {
  matrix_input <- matrix(c(1, 2, 3, 4), nrow = 2)

  expect_equal(pacose:::oppdiag(matrix_input), matrix(c(1, -2, -3, 4), nrow = 2))
})

test_that("Beta2parcor combines paired regression coefficients", {
  beta <- matrix(c(1, 0.8, -0.5,
                   0.5, 1, 0.4,
                   0.2, 0.6, 1), nrow = 3, byrow = TRUE)

  result <- pacose:::Beta2parcor(beta)

  expect_equal(diag(result), rep(1, 3))
  expect_equal(result[1, 2], sqrt(0.8 * 0.5))
  expect_equal(result[1, 3], 0)
})

test_that("isComplete identifies complete and incomplete subgraphs", {
  graph <- igraph::graph_from_literal(1--2, 1--3)

  expect_true(isComplete(graph, c(1, 2)))
  expect_false(isComplete(graph, c(1, 2, 3)))
})

test_that("graph edge counters report intersection and union", {
  first <- igraph::graph_from_literal(1--2, 2--3)
  second <- igraph::graph_from_literal(2--3, 3--4)

  expect_equal(pacose:::nbr(first), 2)
  expect_equal(pacose:::nbri(first, second), 1)
  expect_equal(pacose:::nbru(first, second), 3)
})

test_that("ridge estimation is reproducible with a seed", {
  set.seed(10)
  data <- matrix(rnorm(72), nrow = 24, ncol = 3)
  graph <- igraph::graph_from_literal(1--2, 2--3)

  first <- pacose.ridge(data, graph, k = 3, nlambda = 12, seed = 42)
  second <- pacose.ridge(data, graph, k = 3, nlambda = 12, seed = 42)

  expect_equal(first, second)
})

test_that("INVEST returns a finite matrix and matching constraint error", {
  set.seed(11)
  data <- matrix(rnorm(80), nrow = 20, ncol = 4)
  indices <- matrix(c(0, 1, 1, 2), ncol = 2, byrow = TRUE)

  result <- INVEST_wrapper(data, indices, itermax = 10)

  expect_true(all(is.finite(result$invcovx)))
  expect_equal(
    result$error,
    max(abs(result$invcovx[cbind(indices[, 1] + 1, indices[, 2] + 1)]))
  )
})