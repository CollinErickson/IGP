#' IGP R6 object for fitting hetGP model
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
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {x[1]*sin(2*pi*x[1]) + sqrt(x[1])*exp(x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-2)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_hetGP$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' \donttest{
#' ContourFunctions::cf(function(x) u$predict(x), pts=X1)
#' ContourFunctions::cf(function(x) u$predict(x, se.fit=TRUE)$se, pts=X1)
#' }
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
# \donttest{
#' ContourFunctions::cf(function(x) u$predict(x), pts=rbind(X1, X2))
# }
#' u$delete()
#'
#' n <- 10
#' d <- 1
#' X1 <- runif(n)
#' X1 <- c(X1, X1)
#' f1 <- function(x) {abs(sin(pi*x))+sqrt(x)+rnorm(1,0,.1)}
#' Z1 <- sapply(X1, f1)
#' plot(X1, Z1)
#' h1 <- hetGP::mleHetGP(X=matrix(X1,ncol=1),Z=matrix(Z1, ncol=1), lower=c(.1), upper=c(50))
#' curve(predict(h1, matrix(x, ncol=1))$mean, add=TRUE, col=3)
#' u <- IGP_hetGP$new(X=X1,Z=Z1)
#' plot(X1, Z1, col=2); curve(u$predict(matrix(x,ncol=1)), add=TRUE)
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
IGP_hetGP <- R6::R6Class(
  classname = "IGP_hetGP",
  inherit = IGP_base,
  public = list(
    .init = function(..., kernel="Gaussian", nug.est=self$estimate.nugget, nug=self$nugget0,
                     noiseControl=list(k_theta_g_bounds = c(1, 100),
                                       g_max = 100, g_bounds = c(1e-06, 1)),
                     lower=rep(1e-6, ncol(self$X)), upper=rep(50, ncol(self$X))
                     ) {
      if (!(kernel %in% c("Gaussian", "Matern5_2", "Matern3_2"))) {
        stop('kernel must be in "Gaussian", "Matern5_2", "Matern3_2"')
      }
      if (!nug.est) {noiseControl$g_bounds <- c(nug, nug)}
      # browser()
      m <- hetGP::mleHetGP(X=self$X, Z=matrix(self$Z, ncol=1), covtype=kernel,
                           noiseControl=noiseControl,
                           upper=upper, lower=lower, ...)
      self$mod <- m
    }, #"function to initialize model with data
    .update = function(...) {
      # Can't use update, need Xnew and Znew
      # update(objectself$mod, Xall=self$X, Zall=self$Z, ...)
      # warning("hetGP can't update, will start over")
      self$.delete()
      self$.init()
    }, #"function to add data to model or reestimate params
    .predict = function(XX, se.fit, ...) {
      if (se.fit) {
        preds <- suppressWarnings(predict(object=self$mod, XX))
        list(fit=preds$mean, se.fit=sqrt(preds$sd2))
      } else {
        c(predict(object=self$mod, XX)$mean)
      }
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {sqrt(predict(object=self$mod, XX)$sd2)}, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {predict(object=self$mod, XX)$sd2}, #"function to predict the variance
    .grad = function(XX) NULL, # function to calculate the gradient
    .delete = function(...){self$mod <- NULL}, #"function to delete model beyond simple deletion
    .theta = function() {self$mod$theta}, #"function to get theta, exp(-theta*(x-x)^2)
    .nugget = function() {self$mod$g}, #"function to get nugget
    .s2 = function() {self$mod$nu2_hat},
    .mean = function() {self$mod$beta0} # function that gives mean (constant, other functions not implemented)

  )
)
