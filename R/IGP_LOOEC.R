


#' IGP R6 object for fitting laGP_GauPro model with leave-one-out error correction
#'
#' Class providing object with methods for fitting a GP model.
#' This mixes laGP and GauPro. It fits the model using laGP,
#' then copies the parameters to a GauPro model for prediction.
#' The predicted errors are adjusted by fitting a third GP model
#'  to the leave-one-out absolute t-values.
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
#' u <- IGP_LOOEC_laGP_GauPro$new(X=X1,Z=Z1)
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
#'   \item{Documentation}{For full documentation of each method go to https://github.com/CollinErickson/IGP/}
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL,
#'   estimate.nugget=T, nugget0=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{update(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...)}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
# LOOEC_laGP_GauPro ----
IGP_LOOEC_laGP_GauPro <- R6::R6Class(
  classname = "IGP_LOOEC_laGP_GauPro", inherit = IGP_base,
  public = list(
    package2 = NULL, # Package to fit abstvals model to LOOEC
    # Giving own initialize to get package2
    initialize = function(X=NULL, Z=NULL, package=NULL, corr="gauss", estimate.nugget=TRUE, nugget0=1e-8, package2=NULL, ...) {
      if (!is.null(X)) {self$X <- if (is.matrix(X)) X else matrix(X, ncol=1)} # Add else for 1D data passed as vector
      if (!is.null(Z)) {self$Z <- if (is.matrix(Z)) c(Z) else Z}
      self$package <- package
      self$package2 <- package2
      self$n.at.last.update <- 0
      self$corr <- corr
      self$estimate.nugget <- estimate.nugget
      self$nugget0 <- nugget0
      #if(length(self$X) != 0 & length(self$Z) != 0 & length(self$package) != 0) {
      if(length(self$X) != 0 & length(self$Z) != 0) {
        self$init(...)
      }
    }, # end initialize
    .init = function(..., package2) {
      # Fit model to data with laGP
      self$mod.extra$laGP <- IGP(X=self$X, Z=self$Z, package="laGP",
                                 corr=self$corr,
                                 estimate.nugget=self$estimate.nugget,
                                 nugget0=self$nugget0)
      #self$mod.extra$laGP$init(X=self$X, Z=self$Z, ...)

      # Copy params to GauPro, don't fit, use this for predicting
      self$mod.extra$GauPro <- IGP(X=self$X, Z=self$Z, package="GauPro",
                                   corr=self$corr,
                                   estimate.nugget=FALSE,
                                   nugget0=self$mod.extra$laGP$nugget(),
                                   theta=self$mod.extra$laGP$theta(),
                                   #nug=self$mod.extra$laGP$nug,
                                   param.est=FALSE)

      # Create a third model to model the t values
      # Get warnings for sqrt(-small)
      abstvals <- suppressWarnings({abs(self$mod.extra$GauPro$mod$pred_LOO(se.fit=T)$t)})
      # NaN values are sqrt(-small number), so they should be zero
      # But don't want to make them too small
      abstvals[is.nan(abstvals)] <- 1e-4 #sqrt(.Machine$double.eps)
      # if (!missing(package2)) {self$package2 <- package2}
      if (missing(package2)) {package2 <- "laGP_GauPro"}
      self$package2 <- package2
      self$mod.extra$tmod <- IGP(X=self$X, Z=abstvals, package=self$package2,
                                 corr=self$corr)

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
        self$mod.extra$GauPro$mod$theta <- self$mod.extra$laGP$theta()
        self$mod.extra$GauPro$mod$nug <- self$mod.extra$laGP$nugget()
      }
      self$mod.extra$GauPro$update(Xall=self$X, Zall=self$Z,
                                   no_update=TRUE)

      # Update tmod
      # Get warnings for sqrt(-small number)
      abstvals <- suppressWarnings(abs(self$mod.extra$GauPro$mod$pred_LOO(se.fit=TRUE)$t))
      abstvals[is.nan(abstvals)] <- 1e-4
      self$mod.extra$tmod$update(Xall=self$X, Zall=abstvals, no_update=no_update)

    }, #"function to add data to model or reestimate params
    .predict = function(XX, se.fit, ...) {
      pr <- self$mod.extra$GauPro$.predict(XX=XX, se.fit=se.fit, ...)
      if (se.fit) {
        pt <- self$mod.extra$tmod$predict(XX)
        # Max it up to be at least 1e-4 or the smallest actual t value
        pt <- pmax(1e-4, pt, min(self$mod.extra$tmod$Z))
        pr$se.fit <- pr$se.fit * pt
      }
      pr
    }, #"function to predict at new values
    .predict.se = function(XX, ...) {
      pse <- self$mod.extra$GauPro$.predict.se(XX=XX, ...)
      pt <- self$mod.extra$tmod$predict(XX)
      pt <- pmax(1e-4, pt, min(self$mod.extra$tmod$Z))
      pt * pse
    }, #"function predict the standard error/dev
    .predict.var = function(XX, ...) {
      pv <- self$mod.extra$GauPro$.predict.var(XX=XX, ...)
      pt <- self$mod.extra$tmod$predict(XX)
      pt <- pmax(1e-4, pt, min(self$mod.extra$tmod$Z))
      pt^2 * pv
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
    .mean = NULL # function that gives mean (constant, other functions not implemented)

  )
)
