#' Predict for class IGP
#'
#' @param object Object of class IGP
#' @param XX Points to predict at
#' @param se.fit Whether the standard error prediction should be returned with the mean prediction
#' @param ... Additional parameters
#'
#' @return Prediction from object at XX
#' @export
#'
#' @examples
#' n <- 12
#' x <- matrix(seq(0,1,length.out = n), ncol=1)
#' y <- sin(2*pi*x) + rnorm(n,0,1e-1)
#' gp <- IGP(package='laGP', X=x, Z=y, parallel=FALSE)
#' predict(gp, .448)
predict.IGP <- function(object, XX, se.fit=FALSE, ...) {
  object$predict(XX=XX, se.fit=se.fit)
}
