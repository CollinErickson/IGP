test_that("1D data works", {
  n <- 40
  d <- 2
  n2 <- 20
  f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
  X1 <- matrix(runif(n*d),n,d)
  Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
  gp <- UGP$new(package='laGP',X=X1,Z=Z1)
  expect_that(gp, is_a("UGP"))
  expect_that(gp, is_a("R6"))
  gp$delete()
})
