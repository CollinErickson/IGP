library(R6)
UGP <- R6Class(classname = "UGP",
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
    mod.extra = NULL, #"list", # list to store additional data needed for model
    n.at.last.update = NULL, #"numeric", # count how many in update, must be at end of X
    corr.power = NULL, #"numeric",
    .theta = NULL, #"function",
    .nugget = NULL, #"function",
    estimate.nugget = NULL, #"logical",
    set.nugget = NULL, #"numeric"

    initialize = function(X=NULL, Z=NULL, package=NULL, corr.power=2, ...) {#browser()
      if (!is.null(X)) {self$X <- X}
      if (!is.null(Z)) {self$Z <- if (is.matrix(Z)) c(Z) else Z}
      self$package <- package
      self$n.at.last.update <- 0
      self$corr.power <- corr.power
      #for (item in list(...)) {
      #  self$add(item)
      #}

      if (length(self$package)==0) {
        #message("No package specified Error # 579238572")
      } else if (self$package == "GPfit") {#browser()
        self$.init <- function(...) {
          if (length(self$corr.power) == 0) {
            self$mod <- list(GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=2)))
          } else {
            self$mod <- list(GPfit::GP_fit(self$X, self$Z, corr = list(type="exponential",power=self$corr.power)))
          }
        }
        self$.update <- function(...){
          self$.init()
        }
        self$.predict <- function(XX, se.fit, ...){#browser()
          if (se.fit) {
            preds <- GPfit::predict.GP(self$mod[[1]], XX, se.fit=se.fit)
            list(fit=preds$Y_hat, se.fit=sqrt(preds$MSE))
          } else {
            GPfit::predict.GP(self$mod[[1]], XX)$Y_hat
          }
        }
        self$.predict.se <- function(XX, ...) {sqrt(GPfit::predict.GP(object=self$mod[[1]], xnew=XX, se.fit=T)$MSE)}
        self$.predict.var <- function(XX, ...) {GPfit::predict.GP(object=self$mod[[1]], xnew=XX, se.fit=T)$MSE}
        self$.delete <- function(...){self$mod <- list()}
      } else if (self$package=="laGP") {
        self$.init <- function(...) {#browser()
          da <- laGP::darg(list(mle=TRUE), X=self$X)
          ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=self$Z), silent = T)
          if (inherits(ga.try, "try-error")) {warning("Adding noise to ga in laGP");ga <- laGP::garg(list(mle=TRUE), y=Z+rnorm(length(self$Z),0,1e-2))}
          mod1 <- laGP::newGPsep(X=self$X, Z=self$Z, d=da$start, g=ga$start, dK = TRUE)
          #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
          laGP::jmleGPsep(gpsepi = mod1, drange=c(da$min, da$max),
                                 grange=c(ga$min, ga$max),
                                 dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
          self$mod <- list(mod1)
        }
        self$.update <- function(...) {#browser()
          da <- laGP::darg(list(mle=TRUE), X=self$X)
          ga <- laGP::garg(list(mle=TRUE), y=self$Z)
          n.since.last.update <- nrow(self$X) - self$n.at.last.update
          if (n.since.last.update < 1) {
            message("Can't update, no new X rows")
          } else {
            if (self$n.at.last.update < 10 || n.since.last.update > .25 * self$n.at.last.update) {
              # start over if too many
              self$.delete(...=...)
              self$.init(...=...)
            } else {
              laGP::updateGPsep(gpsepi=self$mod[[1]], X=self$X[-(1:self$n.at.last.update),], Z=self$Z[-(1:self$n.at.last.update)])
            }
          }
          laGP::jmleGPsep(gpsepi = self$mod[[1]], drange=c(da$min, da$max),
                          grange=c(ga$min, ga$max),
                          dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
          }
        self$.predict <- function(XX, se.fit, ...){
          if (se.fit) {
            preds <- laGP::predGPsep(self$mod, XX, lite=TRUE)
            list(fit=preds$mean, se.fit=sqrt(preds$s2))
          } else {
            laGP::predGPsep(self$mod, XX, lite=TRUE)$mean
          }
        }
        self$.predict.se <- function(XX, ...) {sqrt(laGP::predGPsep(self$mod[[1]], XX, lite=TRUE)$s2)}
        self$.predict.var <- function(XX, ...) {laGP::predGPsep(self$mod[[1]], XX, lite=TRUE)$s2}
        self$.delete <- function(...) {
          if (!is.null(self$mod)) {
            laGP::deleteGPsep(self$mod[[1]]);self$mod <- list()
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
          self$mod <- list(mod1)
        }
        self$.update <- function(...) {#browser()
          self$.init(...=...)
        }
        self$.predict <- function(XX, se.fit, ...){#browser()
          capture.output(preds <- with(globalenv(), predict)(self$mod[[1]], XX))
          if (se.fit) {
            list(fit=preds$ZZ.km, se.fit=sqrt(preds$ZZ.ks2))
          } else {
            preds$ZZ.km
          }
        }
        self$.predict.se <- function(XX, ...) {sqrt(with(globalenv(), predict)(self$mod[[1]], XX)$ZZ.ks2)}
        self$.predict.var <- function(XX, ...) {with(globalenv(), predict)(self$mod[[1]], XX)$ZZ.ks2}
        self$.delete <- function(...) {self$mod <- list()}



      } else if (self$package=="mlegp") {
        self$.init <- function(...) {
          co <- capture.output(m <- mlegp::mlegp(X=self$X, Z=self$Z, verbose=0))
          self$mod <- list(m)
        }
        self$.update <- function(...) {
          co <- capture.output(m <- mlegp::mlegp(X=self$X, Z=self$Z, verbose=0))
          self$mod <- list(m)
        }
        self$.predict <- function(XX, se.fit, ...) {
          mlegp::predict.gp(object=self$mod[[1]], newData=XX, se.fit = se.fit)
        }
        self$.predict.se <- function(XX, ...) {mlegp::predict.gp(object=self$mod[[1]], newData=XX, se.fit=T)$se.fit}
        self$.predict.var <- function(XX, ...) {mlegp::predict.gp(object=self$mod[[1]], newData=XX, se.fit=T)$se.fit^2}
        self$.delete <- function(...){self$mod <- list()}


      } else if (self$package=="GauPro") {
        self$.init <- function(...) {
          m <- GauPro::GauPro$new(X=self$X, Z=self$Z, ...)
          self$mod <- list(m)
        }
        self$.update <- function(...) {
          self$mod[[1]]$update(Xall=X, Zall=Z, ...)
        }
        self$.predict <- function(XX, se.fit, ...) {
          if (se.fit) {
            preds <- self$mod[[1]]$pred(XX=XX, se.fit=T)
            list(fit=preds$mean, se.fit=preds$se)
          } else {
            self$mod[[1]]$pred(XX=XX)
          }
        }
        self$.predict.se <- function(XX, ...) {self$mod[[1]]$pred(XX=XX, se.fit=T)$se}
        self$.predict.var <- function(XX, ...) {self$mod[[1]]$pred(XX=XX, se.fit=T)$s2}
        self$.theta <- function() {self$mod[[1]]$theta}
        self$.nugget <- function() {self$mod[[1]]$nug}
        self$.delete <- function(...){self$mod <- list()}


      } else if (self$package=="DiceKriging") {
        self$.init <- function(...) {
          #capture.output(mod1 <- DiceKriging::km(design=X, response=Z, covtype="gauss", nugget.estim=T))
          capture.output(mod1 <- DiceKriging::km(design=self$X, response=self$Z, covtype="gauss", nugget.estim=T))
          self$mod <- list(mod1)
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
              capture.output(DiceKriging::update(object=mod[[1]], newX=self$X[-(1:self$n.at.last.update),], newy=self$Z[-(1:self$n.at.last.update)], nugget.reestim=T))
            } #TRYING TO LET UPDATES BE BIG, ELSE UNCOMMENT THIS PART
          }
        }
        self$.predict <- function(XX, se.fit, ...){
          if (se.fit) {
            preds <- DiceKriging::predict.km(self$mod[[1]], XX, type = "SK", checkNames=F)
            list(fit=preds$mean, se.fit=sqrt(preds$sd))
          } else {
            DiceKriging::predict.km(self$mod[[1]], XX, type = "SK", checkNames=F)$mean
          }
        }
        self$.predict.se <- function(XX, ...) {DiceKriging::predict.km(self$mod[[1]], XX, type = "SK", checkNames=F)$sd}
        self$.predict.var <- function(XX, ...) {(DiceKriging::predict.km(self$mod[[1]], XX, type = "SK", checkNames=F)$sd) ^ 2}
        self$.delete <- function(...) {self$mod <- list()}


      } else if (self$package == "sklearn") {
        #require("rPython")
        self$.init <- function(...) {
          rPython::python.exec('import sys') # These first two lines need to go
          rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
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

          self$mod <- list("GPy model is in Python")
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
          self$mod <- list()
        }
      } else if (self$package == "GPy") {
        #require("rPython")
        self$.init <- function(...) {
          rPython::python.exec('import sys') # These first two lines need to go
          rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
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

          self$mod <- list("GPy model is in Python")
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
          self$mod <- list()
        }
      #} else if (GP.package=='exact') {
      #  predict.GP.SMED <- function(mod,xx) {f(xx)}
      #  init.GP.SMED <- function(X,Y) {}
      #  update.GP.SMED <- function(mod,X,Y) {}
      #  delete.GP.SMED <- function(mod){}
      } else {
        message("Package not recognized Error # 1347344")
      }#;browser()
      if(length(self$X) != 0 & length(self$Z) != 0 & length(self$package) != 0) {
        self$init()
      }
    }, # end initialize
    init = function(X=NULL, Z=NULL, ...) {#browser()
      if (!is.null(X)) {self$X <- X}
      if (!is.null(Z)) {self$Z <- Z}
      if (length(self$X) == 0 | length(self$Z) == 0) {stop("X or Z not set")}
      self$n.at.last.update <- nrow(self$X)
      if (max(self$Z) - min(self$Z) < 1e-8) {warning("Z values are too close, adding noise"); self$Z <- self$Z + rnorm(length(self$Z), 0, 1e-6)}
      #mod <<- list(.init())
      self$.init(...=...)
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...) {#browser()
      if (length(self$n.at.last.update) == 0) {
        init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
      } else {
        if (!is.null(Xall)) {self$X <- Xall} else if (!is.null(Xnew)) {self$X <- rbind(self$X, Xnew)}
        if (!is.null(Zall)) {self$Z <- Zall} else if (!is.null(Znew)) {self$Z <- c(self$Z, Znew)}
        self$.update(...=...)
      }
      self$n.at.last.update <- nrow(self$X)
    }, # end update
    predict = function(XX, se.fit = FALSE, ...) {#browser()
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      self$.predict(XX, se.fit=se.fit, ...=...)
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
    delete = function(...) {
      self$.delete(...=...)
    },
    finalize = function(...) {

    }
  )
)

