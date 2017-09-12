

#' UGP
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
#' \dontrun{
#' n <- 40
#' d <- 2
#' n2 <- 20
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- IGP_GPML$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' }
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
IGP_GPML <- R6::R6Class(classname = "IGP_GPML", inherit = IGP_base,
  public = list(
    matlab_path = NULL,
    .init = function(...) {#browser()

      R.matlab::Matlab$startServer()
      matlab <- R.matlab::Matlab()
      self$mod <- matlab
      isOpen <- open(matlab)
      if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")

      GPML_file_path <- system.file("gpml-matlab-v4.0-2016-10-19", package="UGP")
      # addpath(genpath('C:/Users/cbe117/Documents/R/win-library/3.4/UGP/gpml-matlab-v4.0-2016-10-19'))
      R.matlab::evaluate(matlab, paste0("addpath(genpath('", GPML_file_path, "'));"))
      # set a variable in R and send to MATLAB
      R.matlab::setVariable(matlab, X = self$X)
      R.matlab::setVariable(matlab, Z = self$Z)
      R.matlab::setVariable(matlab, theta = 1)
      R.matlab::evaluate(matlab, 'meanfunc = @meanConst; hyp.mean = [0;];')
      R.matlab::evaluate(matlab, 'covfunc = @covSEard; hyp.cov = [zeros(size(X, 2), 1); 0;]; hyp.lik = log(0.1);')
      R.matlab::evaluate(matlab, 'likfunc = @likGauss')
      # Optimize parameters
      R.matlab::evaluate(matlab, 'hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, X, Z);')
    }, #"function to initialize model with data
    .update = function() { # function to add data to model or reestimate params
      matlab <- self$mod
      R.matlab::setVariable(matlab, X = self$X)
      R.matlab::setVariable(matlab, Z = self$Z)
      R.matlab::evaluate(matlab, "hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, X, Z);")
    },
    .predict = function(XX, se.fit, ...) {#browser()
      R.matlab::setVariable(self$mod, XX = XX)
      R.matlab::evaluate(self$mod, '[mu s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, X, Z, XX);')
      YP <- R.matlab::getVariable(self$mod, 'mu')
      MSEP <- R.matlab::getVariable(self$mod, 's2')
      if (se.fit) {
        cbind(YP$mu, sqrt(MSEP$s2))
      } else {
        YP$mu
      }
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {#browser()
      R.matlab::setVariable(self$mod, XX = XX)
      R.matlab::evaluate(self$mod, '[mu s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, X, Z, XX);')
      YP <- R.matlab::getVariable(self$mod, 'mu')
      MSEP <- R.matlab::getVariable(self$mod, 's2')
      #cbind(YP$YP, sqrt(MSEP$MSEP))
      sqrt(MSEP$s2)
    }, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {#browser()
      R.matlab::setVariable(self$mod, XX = XX)
      R.matlab::evaluate(self$mod, '[mu s2] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, X, Z, XX);')
      # YP <- R.matlab::getVariable(self$mod, 'YP')
      MSEP <- R.matlab::getVariable(self$mod, 's2')
      # cbind(YP$YP, MSEP$MSEP)
      MSEP$s2
    }, #"function to predict the variance
    .grad = NULL, # function to calculate the gradient
    .delete = function() {
      close(self$mod)
      #close(matlab)
    }, #"function to delete model beyond simple deletion
    .theta = function() {
      # R.matlab::getVariable(matlab, 'dmodel.theta')
    }, #"function to get theta, exp(-theta*(x-x)^2)
    .nugget = NULL, #"function to get nugget
    .mean = function() {
      # R.matlab::getVariable(matlab, 'dmodel.beta')
    } # function that gives mean (constant, other functions not implemented)

  )
)



