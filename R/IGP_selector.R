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
IGP <- function(X=NULL, Z=NULL, package=NULL, ...) {#browser()
  if (length(package)==0) {
    stop("No package specified Error # 5792324572")
  } else if (package %in% c("GPfit", "gpfit")) {
    u <- IGP_GPfit$new(X=X, Z=Z, ...)
  } else if (package %in% c("laGP", "laGP")) {
    u <- IGP_laGP$new(X=X, Z=Z, ...)
  } else if (package %in% c("blm","btlm","bcart","bgp","bgpllm","btgp","btgpllm")) {
    u <- IGP_tgp$new(X=X, Z=Z, package=package, ...)
  } else if (package %in% c("mlegp")) {
    u <- IGP_mlegp$new(X=X, Z=Z, ...)
  } else if (package %in% c("GauPro", "gaupro")) {
    u <- IGP_GauPro$new(X=X, Z=Z, ...)
  } else if (package %in% c("GauPro_kernel", "gaupro_kernel")) {
    u <- IGP_GauPro_kernel$new(X=X, Z=Z, ...)
  } else if (tolower(package) %in% c("gaupro_kernel_matern52")) {
    u <- IGP_GauPro_kernel$new(X=X, Z=Z, kernel=GauPro::Matern52, ...)
  } else if (tolower(package) %in% c("gaupro_kernel_matern32")) {
    u <- IGP_GauPro_kernel$new(X=X, Z=Z, kernel=GauPro::Matern32, ...)
  } else if (package %in% c("laGP_GauPro", "lagp_gaupro", "laGP_gaupro", "lagp_GauPro")) {
    u <- IGP_laGP_GauPro$new(X=X, Z=Z, ...)
  } else if (tolower(package) %in% c("lagp_gaupro_kernel")) {
    u <- IGP_laGP_GauPro_kernel$new(X=X, Z=Z, ...)
  } else if (package %in% c("DiceKriging", "dicekriging", "dice", "Dice", "DK", "dk")) {
    u <- IGP_DiceKriging$new(X=X, Z=Z, ...)
  } else if (tolower(package) %in% c("cgp")) {
    u <- IGP_CGP$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "sklearn", "scikit-learn", "scikitlearn")) {
    u <- IGP_sklearn$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "GPy", "gpy")) {
    u <- IGP_GPy$new(X=X, Z=Z, ...)
  } else if (package  %in% c( "DACE", "dace")) {
    # u <- IGP_DACE$new(X=X, Z=Z, ...)
    stop("DACE is currently not available")
  } else if (package  %in% c( "GPML", "gpml")) {
    # u <- IGP_GPML$new(X=X, Z=Z, ...)
    stop("GPML is currently not available")
  } else if (package  %in% c( "LOOEC-laGP_GauPro-laGP")) {
    u <- IGP_LOOEC_laGP_GauPro$new(X=X, Z=Z, package2='laGP', ...)
  } else if (package  %in% c( "LOOEC-laGP_GauPro-GauPro")) {
    u <- IGP_LOOEC_laGP_GauPro$new(X=X, Z=Z, package2='GauPro', ...)
  } else if (tolower(package)  %in% c( "looec-gaupro_kernel")) {
    u <- IGP_LOOEC_GauPro_kernel$new(X=X, Z=Z, package2='GauPro', ...)
  } else {
    stop("Package not recognized")
  }
  u$package <- package
  u
}

