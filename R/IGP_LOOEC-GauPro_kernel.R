
#' IGP R6 object for fitting GauPro model
#'
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' n <- 30
#' d <- 2
#' n2 <- 10
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_LOOEC_GauPro_kernel$new(X=X1,Z=Z1, parallel=FALSE)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' \donttest{
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' }
#' u$delete()
#'
#' \donttest{
#' # 1D example to see difference
#' n <- 9
#' d <- 1
#' n2 <- 20
#' f1 <- function(x) {x^2 * sin(2*pi*x)}
#' X1 <- matrix(seq(0,1,l=n),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-1)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_LOOEC_GauPro_kernel$new(X=X1,Z=Z1, parallel=FALSE)
#' u$plot()
#' u$mod$tmod$plot1D()
#' u$update(Xnew=X2,Znew=Z2)
#' u$plot()
#' u$mod$tmod$plot1D()
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
IGP_LOOEC_GauPro_kernel <- R6::R6Class(
  classname = "IGP_LOOEC_GauPro_kernel",
  inherit = IGP_base,
  public = list(
    .init = function(..., kernel=NULL, theta=NULL) {
      if (!is.null(kernel)) {
        # kernel will be passed in
      } else if (any(c("R6ClassGenerator", "GauPro_kernel")%in% class(self$corr))) {
        kernel <- self$corr
      } else if (self$corr[[1]] == "gauss") {
        kernel <- GauPro::Gaussian$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "matern32") {
        kernel <- GauPro::Matern32$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "matern52") {
        kernel <- GauPro::Matern52$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "exponential") {
        kernel <- GauPro::Exponential$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "periodic") {
        kernel <- GauPro::periodic$new(D=ncol(self$X))
      } else if (self$corr[[1]] == "rationalquadratic") {
        kernel <- GauPro::RatQuad$new(D=ncol(self$X))
      } else {
        stop("Corr/kernel not recognized in IGP_GauPro_kernel")
      }
      if (!is.null(theta)) {kernel$beta <- log(theta, 10)}
      m <- GauPro::GauPro_kernel_model_LOO$new(X=self$X, Z=self$Z, kernel=kernel, nug.est=self$estimate.nugget, nug=self$nugget0, ...)
      self$mod <- m
    }, #"function to initialize model with data
    .update = function(...) {
      self$mod$update(Xall=self$X, Zall=self$Z, ...)
    }, #"function to add data to model or reestimate params
    .predict = function(XX, se.fit, ...) {
      if (se.fit) {
        preds <- self$mod$pred(XX=XX, se.fit=T)
        list(fit=preds$mean, se.fit=preds$se)
      } else {
        c(self$mod$pred(XX=XX))
      }
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$se}, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$s2}, #"function to predict the variance
    .grad = function(XX) {self$mod$grad(XX=XX)}, # function to calculate the gradient
    .delete = function(...){self$mod <- NULL}, #"function to delete model beyond simple deletion
    .theta = function() {10 ^ self$mod$kernel$beta}, #"function to get theta, exp(-theta*(x-x)^2)
    .nugget = function() {self$mod$nug}, #"function to get nugget
    .s2 = function() {self$mod$s2_hat},
    .mean = function() {self$mod$trend$m} # function that gives mean (constant, other functions not implemented)

  )
)
