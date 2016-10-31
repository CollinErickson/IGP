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
#' u <- UGP$new(package='laGP',X=X1,Z=Z1, corr.power=2)
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
UGP <- R6::R6Class(classname = "UGP",
    public = list(
    X = NULL, #"matrix",
    Z = NULL, #"numeric",
    package = NULL, #"character",
    .init = NULL, #"function",
    .update = NULL, #"function",
    .predict = NULL, #"function",
    .predict.se = NULL, #"function",
    .predict.var = NULL, #"function",
    .delete = NULL, #"function",
    mod = NULL, #"list", # First element is model
    mod.extra = list(), #"list", # list to store additional data needed for model
    n.at.last.update = NULL, #"numeric", # count how many in update, must be at end of X
    corr.power = NULL, #"numeric",
    .theta = NULL, #"function",
    .nugget = NULL, #"function",
    estimate.nugget = NULL, #"logical", Should the nugget be estimated?
    set.nugget = NULL, #"numeric" # What value should the nugget be set to? NOT logical
    .mean = NULL, # function that gives mean

    initialize = function(X=NULL, Z=NULL, package=NULL, corr.power=2, estimate.nugget=T, set.nugget=F, ...) {#browser()
      if (!is.null(X)) {self$X <- X}
      if (!is.null(Z)) {self$Z <- if (is.matrix(Z)) c(Z) else Z}
      self$package <- package
      self$n.at.last.update <- 0
      self$corr.power <- corr.power
      self$estimate.nugget <- estimate.nugget
      self$set.nugget <- set.nugget

      if (length(self$package)==0) {
        #message("No package specified Error # 579238572")
      } else if (self$package == "GPfit") {
        self$.init <- function(...) {
          if (!is.null(self$estimate.nugget) || self$set.nugget) {
            warning("GPfit cannot estimate or set the nugget, it picks a stable value")
          }
          if (length(self$corr.power) == 0) {
            self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=2))
          } else {
            self$mod <- GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=self$corr.power))
          }
        }
        self$.update <- function(...){
          self$.init()
        }
        self$.predict <- function(XX, se.fit, ...){#browser()
          if (se.fit) {
            preds <- GPfit::predict.GP(self$mod, XX, se.fit=se.fit)
            list(fit=preds$Y_hat, se.fit=sqrt(preds$MSE))
          } else {
            GPfit::predict.GP(self$mod, XX)$Y_hat
          }
        }
        self$.predict.se <- function(XX, ...) {sqrt(GPfit::predict.GP(object=self$mod, xnew=XX, se.fit=T)$MSE)}
        self$.predict.var <- function(XX, ...) {GPfit::predict.GP(object=self$mod, xnew=XX, se.fit=T)$MSE}
        self$.theta <- function() {10^(self$mod$beta)}
        self$.nugget <- function() {self$mod$delta}
        self$.delete <- function(...){self$mod <- NULL}


      } else if (self$package=="laGP") {
        self$.init <- function(...) {
          da <- laGP::darg(list(mle=TRUE), X=self$X)
          ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z), silent = T)
          if (inherits(ga.try, "try-error")) {
            warning("Adding noise to ga in laGP");
            ga <- laGP::garg(list(mle=TRUE), y=Z+rnorm(length(self$Z),0,1e-2))
          }

          # Follow recommendations for small samples, otherwise use bigger range
          drange <- if (nrow(self$X)<20) c(da$min, da$max) else c(1e-3,1e4) #c(da$min, da$max), # Don't like these small ranges
          grange <- c(ga$min, ga$max)
          mod1 <- laGP::newGPsep(X=self$X, Z=self$Z, d=da$start, g=ga$start, dK = TRUE)
          #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
          mle.out <- laGP::jmleGPsep(gpsepi = mod1,
                                     drange=drange,
                                     grange=grange,
                                     #dab=da$ab, gab=ga$ab, # Will use MLE without these
                                     verb=0, maxit=1000)
          self$mod.extra$theta = as.numeric(1 / mle.out[1,1:ncol(self$X)]) # store theta params
          self$mod.extra$nugget = as.numeric(mle.out[1,ncol(self$X) + 1]) # store nugget
          self$mod <- mod1
        }
        self$.update <- function(...) {
          # Start over if not many points, had problems getting stuck in bad spots early
          if (self$n.at.last.update < 20) {
            self$.delete()
            self$.init(...)
            return()
          }

          da <- laGP::darg(list(mle=TRUE), X=self$X)
          ga <- laGP::garg(list(mle=TRUE), y=self$Z)
          n.since.last.update <- nrow(self$X) - self$n.at.last.update
          if (n.since.last.update < 1) {
            warning("Can't update, no new X rows, but can optimize again")
          } else {
            if (self$n.at.last.update < 10 || n.since.last.update > .25 * self$n.at.last.update) {
              # start over if too many
              self$.delete(...=...)
              self$.init(...=...)
            } else {
              laGP::updateGPsep(gpsepi=self$mod,
                                X=self$X[-(1:self$n.at.last.update),],
                                Z=self$Z[-(1:self$n.at.last.update)])
            }
          }
          drange <- c(1e-3,1e4)
          grange <- c(min(sqrt(.Machine$double.eps),self$mod.extra$nugget), max(1,self$mod.extra$nugget))
          mle.try <- try(mle.out <- laGP::jmleGPsep(gpsepi = self$mod,
                                     #drange=c(da$min, da$max), # Getting rid of these here too
                                     #grange=c(ga$min, ga$max),
                                     drange=drange,
                                     grange=grange, # Had error of nugget starting outside bound
                                     #dab=da$ab, gab=ga$ab,
                                     verb=0, maxit=1000))
          if (inherits(mle.try, "try-error")) {
            # Sometimes gives error: L-BFGS-B needs finite values of 'fn'
            browser()
            warning('Restarting laGP model')
            self$delete()
            self$init(...)
            return()
          }
          # Update stored parameters for when user calls $theta() or $nugget()
          self$mod.extra$theta = as.numeric(1 / mle.out[1,1:ncol(self$X)]) # store theta params
          self$mod.extra$nugget = as.numeric(mle.out[1,ncol(self$X) + 1]) # store nugget
          }
        self$.predict <- function(XX, se.fit, ...){
          if (se.fit) {
            preds <- laGP::predGPsep(self$mod, XX, lite=TRUE)
            list(fit=preds$mean, se.fit=sqrt(preds$s2))
          } else {
            laGP::predGPsep(self$mod, XX, lite=TRUE)$mean
          }
        }
        self$.predict.se <- function(XX, ...) {sqrt(laGP::predGPsep(self$mod, XX, lite=TRUE)$s2)}
        self$.predict.var <- function(XX, ...) {laGP::predGPsep(self$mod, XX, lite=TRUE)$s2}
        self$.theta <- function() {self$mod.extra$theta}
        self$.nugget <- function() {self$mod.extra$nugget}
        self$.delete <- function(...) {
          if (!is.null(self$mod)) {
            laGP::deleteGPsep(self$mod)
            self$mod <- NULL
          }
        }


      } else if (self$package %in% c("blm","btlm","bcart","bgp","bgpllm","btgp","btgpllm")) {
        self$.init <- function(...) {#browser()
          modfunc <-  if (package == "blm") tgp::blm
                      else if (package == "btlm") tgp::btlm
                      else if (package == "bcart") tgp::bcart
                      else if (package == "bgp") tgp::bgp
                      else if (package == "bgpllm") tgp::bgpllm
                      else if (package == "btgp") tgp::btgp
                      else if (package == "btgpllm") tgp::btgpllm
          capture.output(mod1 <- modfunc(self$X, self$Z))
          self$mod <- mod1
        }
        self$.update <- function(...) {#browser()
          self$.init(...=...)
        }
        self$.predict <- function(XX, se.fit, ...){#browser()
          capture.output(preds <- with(globalenv(), predict)(self$mod, XX))
          if (se.fit) {
            list(fit=preds$ZZ.km, se.fit=sqrt(preds$ZZ.ks2))
          } else {
            preds$ZZ.km
          }
        }
        self$.predict.se <- function(XX, ...) {sqrt(with(globalenv(), predict)(self$mod, XX)$ZZ.ks2)}
        self$.predict.var <- function(XX, ...) {with(globalenv(), predict)(self$mod, XX)$ZZ.ks2}
        self$.theta <- function() {rep(NA, ncol(self$X))}
        self$.nugget <- function() {NA}
        self$.delete <- function(...) {self$mod <- NULL}



      } else if (self$package=="mlegp") {
        self$.init <- function(...) {
          temp_nug <- if (is.null(self$estimate.nugget) || self$estimate.nugget == FALSE) NULL
                      else if (self$estimate.nugget == TRUE) 1e-6
          temp_nug_known <- if (is.null(self$set.nugget)) 0 else self$set.nugget
          co <- capture.output(m <- mlegp::mlegp(X=self$X, Z=self$Z, verbose=0,
                                                 nugget = temp_nug,
                                                 nugget.known=temp_nug_known))
          self$mod <- m
        }
        self$.update <- function(...) {
          self$.init(...)
        }
        self$.predict <- function(XX, se.fit, ...) {
          mlegp::predict.gp(object=self$mod, newData=XX, se.fit = se.fit)
        }
        self$.predict.se <- function(XX, ...) {mlegp::predict.gp(object=self$mod, newData=XX, se.fit=T)$se.fit}
        self$.predict.var <- function(XX, ...) {mlegp::predict.gp(object=self$mod, newData=XX, se.fit=T)$se.fit^2}
        self$.theta <- function() {self$mod$beta}
        self$.nugget <- function() {self$mod$nugget}
        self$.delete <- function(...){self$mod <- NULL}


      } else if (self$package=="GauPro") {
        self$.init <- function(...) {
          #m <- GauPro::GauPro$new(X=self$X, Z=self$Z, ...)
          #m <- GauPro::GauPr_Gauss_par$new(X=self$X, Z=self$Z, ...)
          m <- GauPro::GauPro(X=self$X, Z=self$Z, ...)
          self$mod <- m
        }
        self$.update <- function(...) {
          self$mod$update(Xall=self$X, Zall=self$Z, ...)
        }
        self$.predict <- function(XX, se.fit, ...) {
          if (se.fit) {
            preds <- self$mod$pred(XX=XX, se.fit=T)
            list(fit=preds$mean, se.fit=preds$se)
          } else {
            self$mod$pred(XX=XX)
          }
        }
        self$.predict.se <- function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$se}
        self$.predict.var <- function(XX, ...) {self$mod$pred(XX=XX, se.fit=T)$s2}
        self$.theta <- function() {self$mod$theta}
        self$.nugget <- function() {self$mod$nug}
        self$.delete <- function(...){self$mod <- NULL}


      } else if (self$package=="DiceKriging") {
        self$.init <- function(...) {
          #capture.output(mod1 <- DiceKriging::km(design=X, response=Z, covtype="gauss", nugget.estim=T))
          capture.output(mod1 <- DiceKriging::km(design=self$X, response=self$Z, covtype="gauss", nugget.estim=T))
          self$mod <- mod1
        }
        self$.update <- function(...) {#browser()
          n.since.last.update = nrow(self$X) - self$n.at.last.update
          if (n.since.last.update < 1) {
            message("Can't update, no new X rows")
          } else {
            if (self$n.at.last.update < 10 || n.since.last.update > .25 * self$n.at.last.update) {
              # start over if too many
              self$.delete(...=...)
              self$.init(...=...)
            } else {
              capture.output(DiceKriging::update(object=mod, newX=self$X[-(1:self$n.at.last.update),], newy=self$Z[-(1:self$n.at.last.update)], nugget.reestim=T))
            } #TRYING TO LET UPDATES BE BIG, ELSE UNCOMMENT THIS PART
          }
        }
        self$.predict <- function(XX, se.fit, ...){
          if (se.fit) {
            preds <- DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)
            list(fit=preds$mean, se.fit=sqrt(preds$sd))
          } else {
            DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)$mean
          }
        }
        self$.predict.se <- function(XX, ...) {DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)$sd}
        self$.predict.var <- function(XX, ...) {(DiceKriging::predict.km(self$mod, XX, type = "SK", checkNames=F)$sd) ^ 2}
        self$.theta <- function() {self$mod@covariance@range.val}
        self$.nugget <- function() {self$mod@covariance@nugget}
        self$.delete <- function(...) {self$mod <- NULL}


      } else if (self$package == "sklearn") {
        #require("rPython")
        self$.init <- function(...) {
          #rPython::python.exec('import sys') # These first two lines need to go
          #rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
          rPython::python.exec('import numpy as np')
          rPython::python.exec('from sklearn import gaussian_process')
          rPython::python.exec("import warnings")
          rPython::python.exec("warnings.filterwarnings('ignore')")

          rPython::python.assign("inputdim", ncol(self$X))
          rPython::python.assign("X1", (self$X))
          rPython::python.assign("y1", self$Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("gp = gaussian_process.GaussianProcess(                      \
                                theta0=np.asarray([1e-1 for ijk in range(inputdim)]),       \
                                thetaL=np.asarray([1e-4 for ijk in range(inputdim)]),       \
                                thetaU=np.asarray([200 for ijk in range(inputdim)]),        \
                                optimizer='Welch') ")
          rPython::python.exec("gp.fit(X, y)")

          self$mod <- "GPy model is in Python"
        }
        self$.update <- function(...) {
          rPython::python.assign("X1", (self$X))
          rPython::python.assign("y1", self$Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("gp.fit(X, y)")
        }
        self$.predict <- function(XX, se.fit, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(xp, eval_MSE=True)")
          if (se.fit) {
            list(fit=unlist(rPython::python.get("y_pred.tolist()")),
                 se.fit=unlist(rPython::python.get("np.sqrt(sigma2_pred).tolist()")))
          } else {
            unlist(rPython::python.get("y_pred.tolist()"))
          }
        }
        self$.predict.se <- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(xp, eval_MSE=True)")
          unlist(rPython::python.get("np.sqrt(sigma2_pred).tolist()"))
        }
        self$.predict.var <- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(xp, eval_MSE=True)")
          unlist(rPython::python.get("sigma2_pred.tolist()"))
        }
        self$.theta <- function() {rep(NA, ncol(self$X))}
        self$.nugget <- function() {NA}
        self$.delete <- function(...){
          rPython::python.exec('X =  None')
          rPython::python.exec('y =  None')
          rPython::python.exec('xp =  None')
          rPython::python.exec('X1 =  None')
          rPython::python.exec('y1 =  None')
          rPython::python.exec('xp1 =  None')
          rPython::python.exec('y_pred =  None')
          rPython::python.exec('sigma2_pred =  None')
          rPython::python.exec('gp =  None')
          rPython::python.exec('inputdim =  None')
          self$mod <- NULL
        }
      } else if (self$package == "GPy") {
        #require("rPython")
        self$.init <- function(...) {
          #rPython::python.exec('import sys') # These first two lines need to go
          #rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
          rPython::python.exec('import numpy as np')
          rPython::python.exec('import GPy')

          rPython::python.assign("inputdim", ncol(self$X))
          rPython::python.assign("X1", (self$X))
          rPython::python.assign("y1", self$Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("kernel = GPy.kern.RBF(input_dim=inputdim, variance=1., lengthscale=[1. for iii in range(inputdim)],ARD=True)")
          rPython::python.exec("gp = GPy.models.GPRegression(X,y,kernel,normalizer=True)")
          rPython::python.exec("gp.likelihood.variance = 1e-8")
          rPython::python.exec("gp.optimize(messages=False)")
          rPython::python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")

          self$mod <- "GPy model is in Python"
        }
        self$.update <- function(...) {
          rPython::python.assign("X1", (self$X))
          rPython::python.assign("y1", self$Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("gp.set_XY(X = X, Y = y)")
          rPython::python.exec("gp.optimize(messages=False)")
          rPython::python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")
        }
        self$.predict <- function(XX, se.fit, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          if (se.fit) {
            list(fit=unlist(rPython::python.get("y_pred.tolist()")),
                 se.fit=unlist(rPython::python.get("np.sqrt(sigma2_pred).tolist()")))
          } else {
            unlist(rPython::python.get("y_pred.tolist()"))
          }
        }
        self$.predict.se <- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          unlist(rPython::python.get("np.sqrt(sigma2_pred).tolist()"))
        }
        self$.predict.var <- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          unlist(rPython::python.get("sigma2_pred.tolist()"))
        }
        self$.theta <- function() {rep(NA, ncol(self$X))}
        self$.nugget <- function() {NA}
        self$.delete <- function(...){
          rPython::python.exec('X =  None')
          rPython::python.exec('y =  None')
          rPython::python.exec('xp =  None')
          rPython::python.exec('X1 =  None')
          rPython::python.exec('y1 =  None')
          rPython::python.exec('xp1 =  None')
          rPython::python.exec('y_pred =  None')
          rPython::python.exec('sigma2_pred =  None')
          rPython::python.exec('gp =  None')
          rPython::python.exec('kernel =  None')
          rPython::python.exec('inputdim =  None')
          self$mod <- NULL
        }
      #} else if (GP.package=='exact') {
      #  predict.GP.SMED <- function(mod,xx) {f(xx)}
      #  init.GP.SMED <- function(X,Y) {}
      #  update.GP.SMED <- function(mod,X,Y) {}
      #  delete.GP.SMED <- function(mod){}
      } else {
        message("Package not recognized Error # 1347344")
      }
      if(length(self$X) != 0 & length(self$Z) != 0 & length(self$package) != 0) {
        self$init(...)
      }
    }, # end initialize
    init = function(X=NULL, Z=NULL, ...) {#browser()
      if (!is.null(X)) {self$X <- X}
      if (!is.null(Z)) {self$Z <- Z}
      if (length(self$X) == 0 | length(self$Z) == 0) {stop("X or Z not set")}
      self$n.at.last.update <- nrow(self$X)
      if (max(self$Z) - min(self$Z) < 1e-8) {warning("Z values are too close, adding noise"); self$Z <- self$Z + rnorm(length(self$Z), 0, 1e-6)}

      self$.init(...)
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...) {#browser()
      if (self$n.at.last.update == 0) {
        self$init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
      } else {
        if (!is.null(Xall)) {self$X <- Xall} else if (!is.null(Xnew)) {self$X <- rbind(self$X, Xnew)}
        if (!is.null(Zall)) {self$Z <- Zall} else if (!is.null(Znew)) {self$Z <- c(self$Z, Znew)}
        self$.update(...)
      }
      self$n.at.last.update <- nrow(self$X)
    }, # end update
    predict = function(XX, se.fit = FALSE, ...) {#browser()
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
    grad = function (XX) {#browser() # NUMERICAL GRAD IS OVER 10 TIMES SLOWER
      if (!is.matrix(XX)) {
        if (ncol(self$X) == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == ncol(self$X)) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != ncol(self$X)) {stop("Wrong dimension input")}
      }
      grad.func <- function(xx) self$predict(xx)
      grad.apply.func <- function(xx) numDeriv::grad(grad.func, xx)
      grad1 <- apply(XX, 1, grad.apply.func)
      if (ncol(self$X) == 1) return(grad1)
      t(grad1)
    },
    grad_norm = function (XX) {#browser()
      grad1 <- self$grad(XX)
      if (!is.matrix(grad1)) return(abs(grad1))
      apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
    },
    theta = function() {
      self$.theta()
    },
    nugget = function() {
      self$.nugget()
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
    at.max.var = function(X, val=.9) {#browser() #logical if pred var at least 90% of max var
      maxvar = c(self$max.var())
      self$predict.var(X) > val * maxvar
    },
    prop.at.max.var =function(Xlims = matrix(c(0,1), nrow=ncol(self$X), ncol=2, byrow=T), n = 200, val=.9) {#browser()
      maxvar = c(self$max.var())
      X <- apply(Xlims, 1, function(Xlim) {runif(n, Xlim[1], Xlim[2])})
      sum(self$predict.var(X) > val * maxvar) / n
    },
    delete = function(...) {
      self$.delete(...=...)
    },
    finalize = function(...) {
      self$delete() # Mostly for laGP to delete, Python should close connection
    }
  )
)

