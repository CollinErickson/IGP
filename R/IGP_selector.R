#' IGP general function
#'
#' @param package Package to use
#' @param X Design matrix
#' @param Z Response matrix or vector
#' @param ... Arguments passed on to IGP_<package>
#'
#' @return IGP model
#' @export
#'
#' @examples
#' x <- seq(0,1,l=10)
#' y <- abs(sin(2*pi*x))
#' IGP(x,y,'DiceKriging')
IGP <- function(X=NULL, Z=NULL, package=NULL, ...) {
  if (length(package)==0) {
    stop("No package specified Error # 5792324572")
  } else if (package %in% c("GPfit", "gpfit")) {
    u <- IGP_GPfit$new(X=X, Z=Z, ...)
  } else if (package %in% c("laGP", "laGP")) {
    u <- IGP_laGP$new(X=X, Z=Z, ...)
  } else if (package %in% c("blm","btlm","bcart","bgp","bgpllm","btgp","btgpllm")) {
    u <- IGP_tgp$new(X=X, Z=Z, ...)
  } else if (package %in% c("mlegp")) {
    u <- IGP_mlegp$new(X=X, Z=Z, ...)
  } else if (package %in% c("GauPro", "gaupro")) {
    u <- IGP_GauPro$new(X=X, Z=Z, ...)
  } else if (package %in% c("DiceKriging", "dicekriging", "dice", "Dice", "DK", "dk")) {
    u <- IGP_DiceKriging$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "sklearn", "scikit-learn", "scikitlearn")) {
    u <- IGP_sklearn$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "GPy", "gpy")) {
    u <- IGP_GPy$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "DACE", "dace")) {
    u <- IGP_DACE$new(X=X, Z=Z, ...)
  } else {
    stop("Package not recognized")
  }
  u$package <- package
  u
}

