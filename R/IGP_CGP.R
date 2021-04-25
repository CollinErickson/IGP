#' IGP R6 object for fitting CGP model
#'
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @references Ba, S. and V. Roshan Joseph (2012) "Composite Gaussian Process Models for Emulating Expensive Functions". Annals of Applied Statistics, 6, 1838-1860.
#' @examples
#' \donttest{
#' # Takes 17 seconds
#' n <- 20
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_CGP$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' cbind(u$predict(XX1), ZZ1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/IGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, nugget0=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{update(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...)}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_CGP <- R6::R6Class(
  classname = "IGP_CGP",
  inherit = IGP_base,
  public = list(
    .init = function(...) {
      library(CGP) # For some reason gives error if CGP not attached
                   # Can't just use CGP::CGP
      # Suppress warning about recycling array of length 1
      self$mod <- suppressWarnings({CGP::CGP(X=self$X, yobs=self$Z)})
    },
    .update = function(...){
      self$.init()
    }, #"function",
    .predict = function(XX, se.fit, ...){
      if (se.fit) {
        preds <- predict(self$mod, XX)
        list(fit=preds$Yp, se.fit=sqrt(preds$v))
      } else {
        predict(self$mod, XX)$Yp
      }
    }, #"function",
    .predict.se = function(XX, ...) {sqrt(predict(self$mod, XX)$v)}, #"function",
    .predict.var = function(XX, ...) {predict(self$mod, XX)$v}, #"function",
    #.grad = NULL,
    .delete = function(...){self$mod <- NULL}, #"function",
    .theta = function() {self$mod$theta}, #"function",
    .nugget = function() {self$mod$lambda}, #"function",
    .s2 = function() {self$mod$tau2},
    .mean = function() {self$mod$mu} # function that gives mean

  )
)

