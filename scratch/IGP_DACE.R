#' IGP R6 object for fitting DACE model
#'
#' Class providing object with methods for fitting a GP model.
#'
#' DACE was created by Soren N. Lophaven, Hans Bruun Nielsen, Jacob Sondergaard.
#' More information about DACE can be found at http://www2.imm.dtu.dk/projects/dace/.
#'
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
#' u <- IGP_DACE$new(X=X1,Z=Z1)
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
IGP_DACE <- R6::R6Class(classname = "IGP_DACE", inherit = IGP_base,
                        public = list(
                          # matlab_path = "C:\\Users\\cbe117\\School\\DOE\\GP_codes\\DACE\\dace",
                          .init = function(...) {

                            R.matlab::Matlab$startServer()
                            matlab <- R.matlab::Matlab()
                            self$mod <- matlab
                            isOpen <- open(matlab)
                            if (!isOpen) throw("MATLAB server is not running: waited 30 seconds.")

                            # Add DACE folder in IGP to path
                            DACE_file_path <- system.file("dace", package="IGP")
                            R.matlab::evaluate(matlab, paste0("addpath(genpath('", DACE_file_path, "'));"))

                            # set a variable in R and send to MATLAB
                            R.matlab::setVariable(matlab, X = self$X)
                            R.matlab::setVariable(matlab, Z = self$Z)
                            #R.matlab::setVariable(matlab, meanfunc = '')
                            R.matlab::setVariable(matlab, theta = 1)
                            R.matlab::setVariable(matlab, lob = 1e-4)
                            R.matlab::setVariable(matlab, upb = 1e4)
                            R.matlab::evaluate(matlab, 'meanfunc = @regpoly0')
                            R.matlab::evaluate(matlab, 'corrfunc = @corrgauss')
                            # R.matlab::evaluate(matlab, paste0("addpath('",self$matlab_path,"');"))
                            R.matlab::evaluate(matlab, "[dmodel, perf] = dacefit(X, Z, meanfunc, corrfunc, theta, lob, upb);")
                            #R.matlab::evaluate(matlab, "y=20; z=x+y")
                            #z <- R.matlab::getVariable(matlab, "z")
                            #z
                            #close(matlab)

                            #R.matlab::setVariable(matlab, X1=X1)
                            #temp <- R.matlab::evaluate(matlab, "X1 .* X1", capture=TRUE)
                            #temp
                          }, #"function to initialize model with data
                          .update = function() { # function to add data to model or reestimate params
                            matlab <- self$mod
                            R.matlab::setVariable(matlab, X = self$X)
                            R.matlab::setVariable(matlab, Z = self$Z)
                            R.matlab::evaluate(matlab, "[dmodel, perf] = dacefit(X, Z, meanfunc, corrfunc, theta, lob, upb);")
                          },
                          .predict = function(XX, se.fit, ...) {
                            R.matlab::setVariable(self$mod, XX = XX)
                            R.matlab::evaluate(self$mod, '[YP, MSEP] = predictor(XX, dmodel);')
                            YP <- R.matlab::getVariable(self$mod, 'YP')
                            MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                            if (se.fit) {
                              list(fit=YP$YP, se.fit=sqrt(MSEP$MSEP))
                            } else {
                              c(YP$YP)
                            }
                          }, #"function to predict at new values
                          .predict.se = function(XX, ...) {
                            R.matlab::setVariable(self$mod, XX = XX)
                            R.matlab::evaluate(self$mod, '[YP, MSEP] = predictor(XX, dmodel);')
                            # YP <- R.matlab::getVariable(self$mod, 'YP')
                            MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                            #cbind(YP$YP, sqrt(MSEP$MSEP))
                            c(sqrt(MSEP$MSEP))
                          }, #"function predict the standard error/dev
                          .predict.var = function(XX, ...) {
                            R.matlab::setVariable(self$mod, XX = XX)
                            R.matlab::evaluate(self$mod, '[YP, MSEP] = predictor(XX, dmodel);')
                            # YP <- R.matlab::getVariable(self$mod, 'YP')
                            MSEP <- R.matlab::getVariable(self$mod, 'MSEP')
                            # cbind(YP$YP, MSEP$MSEP)
                            c(MSEP$MSEP)
                          }, #"function to predict the variance
                          .grad = NULL, # function to calculate the gradient
                          .delete = function() {
                            close(self$mod)
                            #close(matlab)
                          }, #"function to delete model beyond simple deletion
                          .theta = function() {
                            R.matlab::getVariable(matlab, 'dmodel.theta')
                          }, #"function to get theta, exp(-theta*(x-x)^2)
                          .nugget = NULL, #"function to get nugget
                          .mean = function() {
                            R.matlab::getVariable(matlab, 'dmodel.beta')
                          } # function that gives mean (constant, other functions not implemented)

                        )
)

