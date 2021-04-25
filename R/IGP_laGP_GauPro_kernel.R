
#' IGP R6 object for fitting laGP_GauPro_kernel model
#'
#' Class providing object with methods for fitting a GP model.
#' This mixes laGP and GauPro with a kernel. It fits the model using laGP,
#' then copies the parameters to a GauPro model for prediction.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
# @keywords data, kriging, Gaussian process, regression
#' @return Object of \code{\link{R6Class}} with methods for fitting GP model.
#' @format \code{\link{R6Class}} object.
#' @examples
# \donttest{
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
#' u <- IGP_laGP_GauPro_kernel$new(X=X1,Z=Z1)
#' cbind(u$predict(XX1), ZZ1)
#' u$predict.se(XX1)
#' u$update(Xnew=X2,Znew=Z2)
#' u$predict(XX1)
#' c(u$mod.extra$laGP$theta(), u$mod.extra$laGP$nugget())
#' c(u$mod.extra$GauPro$theta(), u$mod.extra$GauPro$nugget())
#' u$delete()
# }
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
IGP_laGP_GauPro_kernel <- R6::R6Class(
  classname = "IGP_laGP_GauPro_kernel",
  inherit = IGP_base,
  public = list(
    .init = function(...) {
      # Fit model to data with laGP
      self$mod.extra$laGP <- IGP(X=self$X, Z=self$Z, package="laGP",
                                 corr=self$corr,
                                 estimate.nugget=self$estimate.nugget,
                                 nugget0=self$nugget0, ...)
      #self$mod.extra$laGP$init(X=self$X, Z=self$Z, ...)

      # Copy params to GauPro, don't fit, use this for predicting
      kern <- GauPro::Gaussian$new(D=ncol(self$X), beta=log(self$mod.extra$laGP$theta(),10),
                                   s2=self$mod.extra$laGP$s2())
      self$mod.extra$GauPro <- IGP(X=self$X, Z=self$Z, package="GauPro_kernel",
                                   corr=self$corr,
                                   estimate.nugget=FALSE,
                                   nugget0=self$mod.extra$laGP$nugget(),
                                   kernel=kern,
                                   # theta=self$mod.extra$laGP$theta(),
                                   #nug=self$mod.extra$laGP$nug,
                                   param.est=FALSE)
      # laGPs2 <- self$mod.extra$laGP$s2()
      # self$mod.extra$GauPro$mod$kernel$s2 <- laGPs2
      # self$mod.extra$GauPro$mod$kernel$logs2 <- log(laGPs2, 10)
      # self$mod.extra$GauPro$mod$s2_hat <- laGPs2
      #self$mod.extra$GauPro$init(X=self$X, Z=self$Z,
      #                          theta=self$mod.extra$laGP$theta,
      #                           nug=self$mod.extra$laGP$nug)
      self$mod <- "laGP model is mod.extra$laGP, GauPro model is mod.extra$GauPro. This fits with laGP but predicts with GauPro"
    }, #"function to initialize model with data
    .update = function(..., no_update=FALSE) {
      # Update model in laGP
      if (!no_update) { # won't have it update data even not updating params since I don't like when it gives issues, and we pass Xall anyways
        self$mod.extra$laGP$update(Xall=self$X, Zall=self$Z,
                                   no_update=no_update,
                                   ...)
      }
      # Pass GauPro new theta and nugget if it was updated
      if (!no_update) {
        self$mod.extra$GauPro$mod$kernel$beta <- log(self$mod.extra$laGP$theta(), 10)
        self$mod.extra$GauPro$mod$nug <- self$mod.extra$laGP$nugget()
        laGPs2 <- self$mod.extra$laGP$s2()
        self$mod.extra$GauPro$mod$kernel$s2 <- laGPs2
        self$mod.extra$GauPro$mod$kernel$logs2 <- log(laGPs2, 10)
        self$mod.extra$GauPro$mod$s2_hat <- laGPs2
      }
      self$mod.extra$GauPro$update(Xall=self$X, Zall=self$Z,
                                   no_update=TRUE)

    }, #"function to add data to model or reestimate params
    .predict = function(XX, se.fit, ...) {
      self$mod.extra$GauPro$.predict(XX=XX, se.fit=se.fit, ...)
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {
      self$mod.extra$GauPro$.predict.se(XX=XX, ...)
    }, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {
      self$mod.extra$GauPro$.predict.var(XX=XX, ...)
    }, #"function to predict the variance
    .grad = function(XX) {self$mod.extra$GauPro$grad(XX=XX)}, # function to calculate the gradient
    .delete = function(...){
      if (!is.null(self$mod.extra)) {
        self$mod.extra$laGP$delete()
        self$mod.extra$GauPro$delete()
        self$mod.extra <- NULL
      }
      self$mod <- NULL
    }, #"function to delete model beyond simple deletion
    .theta = function() {self$mod.extra$GauPro$theta()}, #"function to get theta, exp(-theta*(x-x)^2)
    .nugget = function() {self$mod.extra$GauPro$nugget()}, #"function to get nugget
    .s2 = function() {self$mod.extra$GauPro$s2()},
    .mean = function() {self$mod.extra$GauPro$trend$m} # function that gives mean (constant, other functions not implemented)

  )
)
