test_that("1D data works", {
  packages = c("laGP", "mlegp", "GauPro", "GPfit", "DiceKriging", "btgp")
  # package = "GauPro"
  for (package in packages) { print(package)
    n <- 40
    d <- 2
    n2 <- 20
    f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
    X1 <- matrix(runif(n*d),n,d)
    Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
    gp <- IGP(package=package,X=X1,Z=Z1)

    # Check predictions work
    nn <- 50
    X2 <- matrix(runif(nn*d),nn,d)
    # Just predict mean
    Z2p <- gp$predict(X2)
    expect_is(Z2p, "numeric", info=package)
    expect_equal(length(Z2p), nn, info=package)
    # Just predict se
    Z2sep <- gp$predict.se(X2)
    expect_is(Z2sep, "numeric", info = package)
    expect_equal(length(Z2sep), nn, info=package)
    # Just predict var
    Z2vp <- gp$predict.var(X2)
    expect_is(Z2vp, "numeric", info=package)
    expect_equal(length(Z2vp), nn, info=package)
    # Just predict mean and se
    Z2p_both <- gp$predict(X2, se.fit = TRUE)
    expect_is(Z2p_both, "list")
    expect_equal(length(Z2p_both$fit), nn, info=package)
    expect_equal(length(Z2p_both$se.fit), nn, info=package)

    # Check classes
    expect_is(gp, "IGP", info=package)
    # expect_that(gp, is_a("IGP_laGP"))
    expect_is(gp, "R6", info=package)

    # Clean-up
    gp$delete()
  }
})
