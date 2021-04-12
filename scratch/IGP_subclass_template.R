#' UGP
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
#' f1 <- function(x) {sin(2*pi*x[1]) + sin(2*pi*x[2])}
#' X1 <- matrix(runif(n*d),n,d)
#' Z1 <- apply(X1,1,f1) + rnorm(n, 0, 1e-3)
#' X2 <- matrix(runif(n2*d),n2,d)
#' Z2 <- apply(X2,1,f1)
#' XX1 <- matrix(runif(10),5,2)
#' ZZ1 <- apply(XX1, 1, f1)
#' u <- UGP2_name123(X=X1,Z=Z1, corr.power=2)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' u$delete()
#' @field X Design matrix
#' @field Z Responses
#' @field N Number of data points
#' @field D Dimension of data
#' @section Methods:
#' \describe{
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/UGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL, corr.power=2,
#'   estimate.nugget=T, set.nugget=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
UGP2_name123 <- R6::R6Class(classname = "UGP2_name123", inherit = UGP2_base,
                          public = list(
                            .init = NULL, #"function to initialize model with data
                            .update = NULL, #"function to add data to model or reestimate params
                            .predict = NULL, #"function to predict at new values
                            .predict.se = NULL, #"function predict the standard error/dev
                            .predict.var = NULL, #"function to predict the variance
                            .grad = NULL, # function to calculate the gradient
                            .delete = NULL, #"function to delete model beyond simple deletion
                            .theta = NULL, #"function to get theta, exp(-theta*(x-x)^2)
                            .nugget = NULL, #"function to get nugget
                            .mean = NULL # function that gives mean (constant, other functions not implemented)

                          )
)
