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
#' u <- IGP(package='laGP',X=X1,Z=Z1, corr="gauss")
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
#'   \item{\code{new(X=NULL, Z=NULL, package=NULL, corr="gauss",
#'   estimate.nugget=T, nugget0=F, ...)}}{This method
#'   is used to create object of this class with \code{X} and \code{Z} as the data.
#'   The package tells it which package to fit the GP model.}
#'   \item{\code{Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...}}{This method
#'   updates the model, adding new data if given, then running optimization again.}}
# IGP_base ----
IGP_base <- R6::R6Class(
  classname = "IGP",
  public = list(
    X = NULL, #"matrix",
    Z = NULL, #"numeric",
    package = NULL, #"character",
    .init = function(...){stop("This function must be overwritten by subclass")}, #"function",
    .update = function(...){stop("This function must be overwritten by subclass")}, #"function",
    .predict = function(...){stop("This function must be overwritten by subclass")}, #"function",
    .predict.se = function(...){stop("This function must be overwritten by subclass")}, #"function",
    .predict.var = function(...){stop("This function must be overwritten by subclass")}, #"function",
    #.grad = function(...){stop("This function must be overwritten by subclass")},
    .delete = function(...){self$mod <- NULL}, #"function",
    #.theta = function(...){stop("This function must be overwritten by subclass")}, #"function",
    #.nugget = function(...){stop("This function must be overwritten by subclass")}, #"function",
    #.mean = function(...){stop("This function must be overwritten by subclass")}, # function that gives mean
    mod = NULL, #"list", # First element is model
    mod.extra = list(), #"list", # list to store additional data needed for model
    n.at.last.update = NULL, #"numeric", # count how many in update, must be at end of X
    corr = NULL, #"numeric",
    estimate.nugget = NULL, #"logical", Should the nugget be estimated?
    nugget0 = NULL, #"numeric" # What value should the nugget be set to? NOT logical.
    # If estimate.nugget==TRUE, then it's the starting value

    initialize = function(X=NULL, Z=NULL, package=NULL, corr="gauss", estimate.nugget=TRUE, nugget0=1e-8, ...) {
      if (!is.null(X)) {self$X <- if (is.matrix(X)) X else matrix(X, ncol=1)} # Add else for 1D data passed as vector
      if (!is.null(Z)) {self$Z <- if (is.matrix(Z)) c(Z) else Z}
      self$package <- package
      self$n.at.last.update <- 0
      self$corr <- corr
      self$estimate.nugget <- estimate.nugget
      self$nugget0 <- nugget0

      #if(length(self$X) != 0 & length(self$Z) != 0 & length(self$package) != 0) {
      if(length(self$X) != 0 & length(self$Z) != 0) {
        self$init(...)
      }
    }, # end initialize
    init = function(X=NULL, Z=NULL, ...) {
      if (!is.null(X)) {self$X <- X}
      if (!is.null(Z)) {self$Z <- Z}
      if (length(self$X) == 0 | length(self$Z) == 0) {stop("X or Z not set")}
      self$n.at.last.update <- nrow(self$X)
      if (max(self$Z) - min(self$Z) < 1e-8) {
        warning("Z values are too close, adding noise")
        self$Z <- self$Z + rnorm(length(self$Z), 0, 1e-6)
      }

      self$.init(...)
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...) {
      if (self$n.at.last.update == 0) {
        #self$init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
        x <- if(!is.null(Xall)) Xall else Xnew
        z <- if (!is.null(Zall)) Zall else Znew
        self$init(X = x, Z = z)
      } else {
        if (!is.null(Xall)) {self$X <- Xall} else if (!is.null(Xnew)) {self$X <- rbind(self$X, Xnew)}
        if (!is.null(Zall)) {self$Z <- Zall} else if (!is.null(Znew)) {self$Z <- c(self$Z, Znew)}
        self$.update(...)
      }
      self$n.at.last.update <- nrow(self$X)
    }, # end update
    predict = function(XX, se.fit = FALSE, ...) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      self$.predict(XX, se.fit=se.fit, ...)
    },
    predict.se = function(XX, ...) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      self$.predict.se(XX, ...=...)
    },
    predict.var = function(XX, ...) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      self$.predict.var(XX, ...=...)
    },
    grad = function (XX, num=FALSE) { # NUMERICAL GRAD IS OVER 10 TIMES SLOWER
      if (!is.matrix(XX)) {
        if (ncol(self$X) == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == ncol(self$X)) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != ncol(self$X)) {stop("Wrong dimension input")}
      }
      if (is.null(self$.grad) | num) { # if no method, use numerical
        #print('using num')
        self$grad_num(XX)
      } else {#print('using package')
        self$.grad(XX)
      }
    },
    grad_num = function (XX) {
      grad.func <- function(xx) self$predict(xx)
      grad.apply.func <- function(xx) numDeriv::grad(grad.func, xx)
      grad1 <- apply(XX, 1, grad.apply.func)
      if (ncol(self$X) == 1) return(grad1)
      t(grad1)
    },
    grad_from_theta = function(XX, theta) {
      if (missing(theta)) {
        theta <- self$theta()
        if (is.null(theta)) {
          stop("Need theta for grad_from_theta")
        }
      }
      mu <- self$mean()
      D <- ncol(self$X)
      N <- nrow(self$X)
      if (!is.matrix(XX)) {
        if (D == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == D) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != D) {stop("Wrong dimension input")}
      }
      # kx.xx <- self$corr_func(self$X, XX, theta=self$theta)
      kx.xx <- GauPro::corr_gauss_matrix(self$X, XX, theta)
      Kx <- GauPro::corr_gauss_matrix_symC(self$X, theta)
      Kx_nug <- Kx + diag(self$nugget(), nrow(Kx))
      Kinv_Z_minus_mu <- solve(Kx_nug, self$Z - mu)

      grad1 <-   vapply(1:nrow(XX),
                        Vectorize(
                          function(k) {
                            t(-2 * outer(1:N, 1:D,
                                         Vectorize(function(i,j) {
                                           theta[j] * (XX[k, j] - self$X[i, j]) * kx.xx[i, k]
                                         }))
                            )  %*%Kinv_Z_minus_mu
                          }
                        )
                        , numeric(D)
      )
      if (D == 1) return(grad1)
      t(grad1)
    },
    grad_norm = function (XX) {
      grad1 <- self$grad(XX)
      if (!is.matrix(grad1)) return(abs(grad1))
      apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
    },
    sample = function(XX, n=1) {
      if (length(XX) != ncol(self$X)) {stop("Can only sample one point at a time right now error 23537898")}
      XX.pred <- self$predict(XX=XX, se.fit=T)
      rnorm(n=n, mean=XX.pred$fit, sd=XX.pred$se.fit)
    },
    theta = function() {
      self$.theta()
    },
    nugget = function() {
      self$.nugget()
    },
    s2 = function() {
      self$.s2()
    },
    mean = function() {
      if (!is.null(self$.mean)) {
        self$.mean()
      } else {
        self$predict(matrix(rep(max(abs(self$X)) * 10,ncol(self$X)), nrow=1))
      }
    },
    max.var = function() {
      self$predict.var(matrix(rep(max(abs(self$X)) * 10,ncol(self$X)), nrow=1))
    },
    at.max.var = function(X, val=.9) { #logical if pred var at least 90% of max var
      maxvar = c(self$max.var())
      self$predict.var(X) > val * maxvar
    },
    prop.at.max.var =function(Xlims = matrix(c(0,1), nrow=ncol(self$X), ncol=2, byrow=T), n = 200, val=.9) {
      maxvar = c(self$max.var())
      X <- apply(Xlims, 1, function(Xlim) {runif(n, Xlim[1], Xlim[2])})
      sum(self$predict.var(X) > val * maxvar) / n
    },
    plot = function() {
      minx <- min(self$X)
      maxx <- max(self$X)
      minxeval <- minx - .03 * (maxx - minx)
      maxxeval <- maxx + .03 * (maxx - minx)
      if (ncol(self$X) == 1) {
        XX <- matrix(seq(minxeval,maxxeval,length.out = 300), ncol=1)
        pp <- self$predict(XX=XX, se.fit=TRUE)
        pm <- pp$fit
        ps <- pp$se.fit
        phigh <- pm + 2 * ps
        plow  <- pm - 2 * ps
        plot(XX, pm, col='white', ylim=c(min(plow), max(phigh)),
             xlab="X", ylab="Z")
        points(XX, phigh, type='l', col=2, lwd=2)
        points(XX, plow, type='l', col=2, lwd=2)
        points(XX, pm, type='l', lwd=3)
        points(self$X, self$Z, pch=19, cex=2)
      } else {
        stop("No plot method for higher than 1D")
      }
    },
    delete = function(...) {
      self$.delete(...=...)
    },
    finalize = function(...) {
      self$delete() # Mostly for laGP to delete, Python should close connection
    }
  )
)

