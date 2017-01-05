UGP2 <- function(..., package=NULL) {

  if (length(package)==0) {
    stop("No package specified Error # 5792324572")
  } else if (self$package == "GPfit") {
    u <- UGP2_GPfit$new(...)
  } else if (package=="laGP") {
    u <- UGP2_laGP$new(...)
  } else if (package %in% c("blm","btlm","bcart","bgp","bgpllm","btgp","btgpllm")) {
    u <- UGP2_tgp$new(...)
  } else if (package=="mlegp") {
    u <- UGP2_mlegp$new(...)
  } else if (package=="GauPro") {
    u <- UGP2_GauPro$new(...)
  } else if (package=="DiceKriging") {
    u <- UGP2_DiceKriging$new(...)
  } else if (package == "sklearn") {
    u <- UGP2_sklearn$new(...)
  } else if (package == "GPy") {
    u <- UGP2_GPy$new(...)
  } else {
    stop("Package not recognized")
  }
  u$package <- package
  u
}

