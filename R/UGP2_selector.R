#' UGP2 general function
#'
#' @param X=X, Z=Z, ... Passed on to specific package model
#' @param package Package to use
#' @param X Design matrix
#' @param Z Response matrix or vector
#'
#' @return UGP2 model
#' @export
#'
#' @examples
#' x <- seq(0,1,l=10)
#' y <- abs(sin(2*pi*x))
#' UGP2(x,y,'DiceKriging')
UGP2 <- function(X=NULL, Z=NULL, package=NULL, ...) {
  if (length(package)==0) {
    stop("No package specified Error # 5792324572")
  } else if (package %in% c("GPfit", "gpfit")) {
    u <- UGP2_GPfit$new(X=X, Z=Z, ...)
  } else if (package %in% c("laGP", "laGP")) {
    u <- UGP2_laGP$new(X=X, Z=Z, ...)
  } else if (package %in% c("blm","btlm","bcart","bgp","bgpllm","btgp","btgpllm")) {
    u <- UGP2_tgp$new(X=X, Z=Z, ...)
  } else if (package %in% c("mlegp")) {
    u <- UGP2_mlegp$new(X=X, Z=Z, ...)
  } else if (package %in% c("GauPro", "gaupro")) {
    u <- UGP2_GauPro$new(X=X, Z=Z, ...)
  } else if (package %in% c("DiceKriging", "dicekriging", "dice", "Dice", "DK", "dk")) {
    u <- UGP2_DiceKriging$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "sklearn", "scikit-learn", "scikitlearn")) {
    u <- UGP2_sklearn$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "GPy", "gpy")) {
    u <- UGP2_GPy$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "DACE", "dace")) {
    u <- UGP2_DACE$new(X=X, Z=Z, ...)
  } else {
    stop("Package not recognized")
  }
  u$package <- package
  u
}

