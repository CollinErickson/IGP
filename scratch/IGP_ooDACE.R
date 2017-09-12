#' UGP ooDACE model
#'
#' Class providing object with methods for fitting a GP model
#'
#' @docType class
#' @importFrom R6 R6Class
# @export
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
#' u <- IGP_ooDACE$new(X=X1,Z=Z1)
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
IGP_ooDACE <- R6::R6Class(classname = "IGP_ooDACE", inherit = IGP_base,
  public = list(
    .init = function(...) {browser()

      R.matlab::Matlab$startServer()
      matlab <- R.matlab::Matlab()
      self$mod <- matlab
      isOpen <- open(matlab)
      if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")

      # set a variable in R and send to MATLAB
      R.matlab::setVariable(matlab, X = self$X)
      R.matlab::setVariable(matlab, inputdim = ncol(self$X))
      R.matlab::setVariable(matlab, Z = self$Z)
      #R.matlab::setVariable(matlab, meanfunc = '')
      R.matlab::setVariable(matlab, theta = 1)
      R.matlab::setVariable(matlab, lob = 1e-4)
      R.matlab::setVariable(matlab, upb = 1e4)
      R.matlab::evaluate(matlab, 'meanfunc = @regpoly0')
      R.matlab::evaluate(matlab, 'corrfunc = @corrgauss')
      R.matlab::evaluate(matlab, "addpath(genpath('C:\\Users\\cbe117\\School\\DOE\\GP_codes\\ooDACE'));")
      # R.matlab::evaluate(matlab, "addpath('C:\\Users\\cbe117\\School\\DOE\\GP_codes\\ooDACE');")
      # R.matlab::evaluate(matlab, "addpath('C:\\Users\\cbe117\\School\\DOE\\GP_codes\\ooDACE\\optimizers\');")
      #R.matlab::evaluate(matlab, "[dmodel, perf] = DACEfit(X, Z, meanfunc, corrfunc, theta, lob, upb);")
      browser()
      R.matlab::evaluate(matlab, "opts = Kriging.getDefaultOptions();")
      R.matlab::evaluate(matlab, "opts.hpBounds = [lob ; upb]; % hyperparameter optimization bounds")
      R.matlab::evaluate(matlab, "opts.hpOptimizer = SQPLabOptimizer( inputdim, 1 );")
      if (self$estimate.nugget) {
        R.matlab::evaluate(matlab, "opts.lambda0 = 0;")
        R.matlab::evaluate(matlab, "opts.lambdaBounds = [-5; 5];")
      }
      R.matlab::evaluate(matlab, "")
      R.matlab::evaluate(matlab, "")
      R.matlab::evaluate(matlab, "")
      R.matlab::evaluate(matlab, "k = Kriging( opts, theta, meanfunc, corrfunc);")
      R.matlab::evaluate(matlab, "k = k.fit( X, Z);")
      #R.matlab::evaluate(matlab, "y=20; z=x+y")
      #z <- R.matlab::getVariable(matlab, "z")
      #z
      #close(matlab)

      #R.matlab::setVariable(matlab, X1=X1)
      #temp <- R.matlab::evaluate(matlab, "X1 .* X1", capture=TRUE)
      #temp
    }, #"function to initialize model with data
    .update = function (...) { #function to add data to model or reestimate params
      matlab <- self$mod
      R.matlab::setVariable(matlab, X = self$X)
      R.matlab::setVariable(matlab, Z = self$Z)
      R.matlab::evaluate(matlab, "k = k.fit( X, Z);")
    },
    .predict = function(XX, se.fit, ...) {browser()
      R.matlab::setVariable(self$mod, XX = XX)
      R.matlab::evaluate(self$mod, '[YP, MSEP] = k.predict(XX);')
      YP <- R.matlab::getVariable(self$mod, 'YP')
      MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
      if (se.fit) {
        cbind(YP$YP, sqrt(MSEP$MSEP))
      } else {
        YP$YP
      }
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {browser()
      R.matlab::setVariable(self$mod, XX = XX)
      R.matlab::evaluate(self$mod, '[YP, MSEP] = k.predict(XX);')
      YP <- R.matlab::getVariable(self$mod, 'YP')
      MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
      YP$YP
      cbind(YP$YP, sqrt(MSEP$MSEP))
    }, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {browser()
      R.matlab::setVariable(self$mod, XX = XX)
      R.matlab::evaluate(self$mod, '[YP, MSEP] = k.predict(XX);')
      YP <- R.matlab::getVariable(self$mod, 'YP')
      MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
      cbind(YP$YP, MSEP$MSEP)
    }, #"function to predict the variance
    .grad = NULL, # function to calculate the gradient
    .delete = function() {
      close(self$mod)
      #close(matlab)
    }, #"function to delete model beyond simple deletion
    .theta = function() {
      R.matlab::getVariable(matlab, 'k.getHyperparameters')
    }, #"function to get theta, exp(-theta*(x-x)^2)
    .nugget = NULL, #"function to get nugget
    .mean = function() {
      R.matlab::getVariable(matlab, 'methods(k)') # can't do this yet
    } # function that gives mean (constant, other functions not implemented)

  )
)
