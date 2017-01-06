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
  } else if (package == "GPfit") {
    u <- UGP2_GPfit$new(X=X, Z=Z, ...)
  } else if (package=="laGP") {
    u <- UGP2_laGP$new(X=X, Z=Z, ...)
  } else if (package %in% c("blm","btlm","bcart","bgp","bgpllm","btgp","btgpllm")) {
    u <- UGP2_tgp$new(X=X, Z=Z, ...)
  } else if (package=="mlegp") {
    u <- UGP2_mlegp$new(X=X, Z=Z, ...)
  } else if (package=="GauPro") {
    u <- UGP2_GauPro$new(X=X, Z=Z, ...)
  } else if (package=="DiceKriging") {
    u <- UGP2_DiceKriging$new(X=X, Z=Z, ...)
  } else if (package == "sklearn") {
    u <- UGP2_sklearn$new(X=X, Z=Z, ...)
  } else if (package == "GPy") {
    u <- UGP2_GPy$new(X=X, Z=Z, ...)
  } else {
    stop("Package not recognized")
  }
  u$package <- package
  u
}

