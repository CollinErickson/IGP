UGP <- setRefClass("UGP",
  fields = list(
    X = "matrix",
    Z = "numeric",
    package = "character",
    .init = "function",
    .update = "function",
    .predict = "function",
    .predict.se = "function",
    .predict.var = "function",
    .delete = "function",
    mod = "list", # First element is model
    mod.extra = "list", # list to store additional data needed for model
    n.at.last.update = "numeric", # count how many in update, must be at end of X
    corr.power = "numeric",
    .theta = "function",
    .nugget = "function",
    estimate.nugget = "logical",
    set.nugget = "numeric"
  ),
  methods = list(
    initialize = function(X=NULL, Z=NULL, ...) {#browser()
      if (!is.null(X)) {X <<- X}
      if (!is.null(Z)) {Z <<- if (is.matrix(Z)) c(Z) else Z}
      callSuper(...)

      if (length(package)==0) {
        #message("No package specified Error # 579238572")
      } else if (package == "GPfit") {
        .init <<- function(...) {
          if (length(corr.power) == 0) {
            mod <<- list(GPfit::GP_fit(X, Z, corr = list(type="exponential",power=2)))
          } else {
            mod <<- list(GPfit::GP_fit(X, Z, corr = list(type="exponential",power=corr.power)))
          }
        }
        .update <<- function(...){
          .init()
        }
        .predict <<- function(XX, se.fit, ...){
          if (se.fit) {
            preds <- GPfit::predict.GP(mod[[1]], XX, se.fit=se.fit)
            list(fit=preds$Y_hat, se.fit=sqrt(preds$MSE))
          } else {
            GPfit::predict.GP(mod[[1]], XX)$Y_hat
          }
        }
        .predict.se <<- function(XX, ...) {sqrt(GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE)}
        .predict.var <<- function(XX, ...) {GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE}
        .delete <<- function(...){mod <<- list()}
      } else if (package=="laGP") {
        .init <<- function(...) {#browser()
          da <- laGP::darg(list(mle=TRUE), X=X)
          ga.try <- try(ga <- laGP::garg(list(mle=TRUE), y=Z), silent = T)
          if (inherits(ga.try, "try-error")) {warning("Adding noise to ga in laGP");ga <- laGP::garg(list(mle=TRUE), y=Z+rnorm(length(Z),0,1e-2))}
          mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=ga$start, dK = TRUE)
          #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
          laGP::jmleGPsep(gpsepi = mod1, drange=c(da$min, da$max),
                                 grange=c(ga$min, ga$max),
                                 dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
          mod <<- list(mod1)
        }
        .update <<- function(...) {#browser()
          da <- laGP::darg(list(mle=TRUE), X=X)
          ga <- laGP::garg(list(mle=TRUE), y=Z)
          n.since.last.update = nrow(X) - n.at.last.update
          if (n.since.last.update < 1) {
            message("Can't update, no new X rows")
          } else {
            if (n.at.last.update < 10 || n.since.last.update > .25 * n.at.last.update) {
              # start over if too many
              .delete(...=...)
              .init(...=...)
            } else {
              laGP::updateGPsep(gpsepi=mod[[1]], X=X[-(1:n.at.last.update),], Z=Z[-(1:n.at.last.update)])
            }
          }
          laGP::jmleGPsep(gpsepi = mod[[1]], drange=c(da$min, da$max),
                          grange=c(ga$min, ga$max),
                          dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
          }
        .predict <<- function(XX, se.fit, ...){
          if (se.fit) {
            preds <- laGP::predGPsep(mod, XX, lite=TRUE)
            list(fit=preds$mean, se.fit=sqrt(preds$s2))
          } else {
            laGP::predGPsep(mod, XX, lite=TRUE)$mean
          }
        }
        .predict.se <<- function(XX, ...) {sqrt(laGP::predGPsep(mod, XX, lite=TRUE)$s2)}
        .predict.var <<- function(XX, ...) {laGP::predGPsep(mod, XX, lite=TRUE)$s2}
        .delete <<- function(...) {laGP::deleteGPsep(mod[[1]]);mod <<- list()}


      } else if (package %in% c("blm","btlm","bcart","bgp","bgpllm","btgp","btgpllm")) {
        .init <<- function(...) {#browser()
          modfunc <-  if (package == "blm") tgp::blm
                      else if (package == "btlm") tgp::btlm
                      else if (package == "bcart") tgp::bcart
                      else if (package == "bgp") tgp::bgp
                      else if (package == "bgpllm") tgp::bgpllm
                      else if (package == "btgp") tgp::btgp
                      else if (package == "btgpllm") tgp::btgpllm
          capture.output(mod1 <- modfunc(X, Z))
          mod <<- list(mod1)
        }
        .update <<- function(...) {#browser()
          .init(...=...)
        }
        .predict <<- function(XX, se.fit, ...){#browser()
          capture.output(preds <- with(globalenv(), predict)(mod[[1]], XX))
          if (se.fit) {
            list(fit=preds$ZZ.km, se.fit=sqrt(preds$ZZ.ks2))
          } else {
            preds$ZZ.km
          }
        }
        .predict.se <<- function(XX, ...) {sqrt(with(globalenv(), predict)(mod[[1]], XX)$ZZ.ks2)}
        .predict.var <<- function(XX, ...) {with(globalenv(), predict)(mod[[1]], XX)$ZZ.ks2}
        .delete <<- function(...) {mod <<- list()}



      } else if (package=="mlegp") {
        .init <<- function(...) {
          co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
          mod <<- list(m)
        }
        .update <<- function(...) {
          co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
          mod <<- list(m)
        }
        .predict <<- function(XX, se.fit, ...) {
          mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit = se.fit)
        }
        .predict.se <<- function(XX, ...) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit}
        .predict.var <<- function(XX, ...) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit^2}
        .delete <<- function(...){mod <<- list()}


      } else if (package=="GauPro") {
        .init <<- function(...) {
          m <- GauPro::GauPro$new(X=X, Z=Z, ...)
          mod <<- list(m)
        }
        .update <<- function(...) {
          mod[[1]]$update(Xall=X, Zall=Z, ...)
        }
        .predict <<- function(XX, se.fit, ...) {
          if (se.fit) {
            preds <- mod[[1]]$pred(XX=XX, se.fit=T)
            list(fit=preds$mean, se.fit=preds$se)
          } else {
            mod[[1]]$pred(XX=XX)
          }
        }
        .predict.se <<- function(XX, ...) {mod[[1]]$pred(XX=XX, se.fit=T)$se}
        .predict.var <<- function(XX, ...) {mod[[1]]$pred(XX=XX, se.fit=T)$s2}
        .theta <<- function() {mod[[1]]$theta}
        .nugget <<- function() {mod[[1]]$nug}
        .delete <<- function(...){mod <<- list()}


      } else if (package=="DiceKriging") {
        .init <<- function(...) {
          #capture.output(mod1 <- DiceKriging::km(design=X, response=Z, covtype="gauss", nugget.estim=T))
          capture.output(mod1 <- DiceKriging::km(design=X, response=Z, covtype="gauss", nugget.estim=T))
          mod <<- list(mod1)
        }
        .update <<- function(...) {#browser()
          n.since.last.update = nrow(X) - n.at.last.update
          if (n.since.last.update < 1) {
            message("Can't update, no new X rows")
          } else {
            if (n.at.last.update < 10 || n.since.last.update > .25 * n.at.last.update) {
              # start over if too many
              .delete(...=...)
              .init(...=...)
            } else {
              capture.output(DiceKriging::update(object=mod[[1]], newX=X[-(1:n.at.last.update),], newy=Z[-(1:n.at.last.update)], nugget.reestim=T))
            } #TRYING TO LET UPDATES BE BIG, ELSE UNCOMMENT THIS PART
          }
        }
        .predict <<- function(XX, se.fit, ...){
          if (se.fit) {
            preds <- DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)
            list(fit=preds$mean, se.fit=sqrt(preds$sd))
          } else {
            DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)$mean
          }
        }
        .predict.se <<- function(XX, ...) {DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)$sd}
        .predict.var <<- function(XX, ...) {(DiceKriging::predict.km(mod[[1]], XX, type = "SK", checkNames=F)$sd) ^ 2}
        .delete <<- function(...) {mod <<- list()}


      } else if (package == "sklearn") {
        #require("rPython")
        .init <<- function(...) {
          rPython::python.exec('import sys') # These first two lines need to go
          rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
          rPython::python.exec('import numpy as np')
          rPython::python.exec('from sklearn import gaussian_process')
          rPython::python.exec("import warnings")
          rPython::python.exec("warnings.filterwarnings('ignore')")

          rPython::python.assign("inputdim", ncol(X))
          rPython::python.assign("X1", (X))
          rPython::python.assign("y1", Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("gp = gaussian_process.GaussianProcess(                      \
                                theta0=np.asarray([1e-1 for ijk in range(inputdim)]),       \
                                thetaL=np.asarray([1e-4 for ijk in range(inputdim)]),       \
                                thetaU=np.asarray([200 for ijk in range(inputdim)]),        \
                                optimizer='Welch') ")
          rPython::python.exec("gp.fit(X, y)")

          mod <<- list("GPy model is in Python")
        }
        .update <<- function(...) {
          rPython::python.assign("X1", (X))
          rPython::python.assign("y1", Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("gp.fit(X, y)")
        }
        .predict <<- function(XX, se.fit, ...) {
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
        .predict.se <<- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(xp, eval_MSE=True)")
          unlist(rPython::python.get("np.sqrt(sigma2_pred).tolist()"))
        }
        .predict.var <<- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(xp, eval_MSE=True)")
          unlist(rPython::python.get("sigma2_pred.tolist()"))
        }
        .delete <<- function(...){
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
          mod <<- list()
        }
      } else if (package == "GPy") {
        #require("rPython")
        .init <<- function(...) {
          rPython::python.exec('import sys') # These first two lines need to go
          rPython::python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
          rPython::python.exec('import numpy as np')
          rPython::python.exec('import GPy')

          rPython::python.assign("inputdim", ncol(X))
          rPython::python.assign("X1", (X))
          rPython::python.assign("y1", Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("kernel = GPy.kern.RBF(input_dim=inputdim, variance=1., lengthscale=[1. for iii in range(inputdim)],ARD=True)")
          rPython::python.exec("gp = GPy.models.GPRegression(X,y,kernel,normalizer=True)")
          rPython::python.exec("gp.likelihood.variance = 1e-8")
          rPython::python.exec("gp.optimize(messages=False)")
          rPython::python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")

          mod <<- list("GPy model is in Python")
        }
        .update <<- function(...) {
          rPython::python.assign("X1", (X))
          rPython::python.assign("y1", Z)
          rPython::python.exec('X =  np.matrix(X1)')
          rPython::python.exec('y = np.matrix(y1).reshape((-1,1))')
          rPython::python.exec("gp.set_XY(X = X, Y = y)")
          rPython::python.exec("gp.optimize(messages=False)")
          rPython::python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")
        }
        .predict <<- function(XX, se.fit, ...) {
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
        .predict.se <<- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          unlist(rPython::python.get("np.sqrt(sigma2_pred).tolist()"))
        }
        .predict.var <<- function(XX, ...) {
          rPython::python.assign("xp1", XX)
          rPython::python.exec("xp = np.asmatrix(xp1)")
          rPython::python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          unlist(rPython::python.get("sigma2_pred.tolist()"))
        }
        .delete <<- function(...){
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
          mod <<- list()
        }
      #} else if (GP.package=='exact') {
      #  predict.GP.SMED <- function(mod,xx) {f(xx)}
      #  init.GP.SMED <- function(X,Y) {}
      #  update.GP.SMED <- function(mod,X,Y) {}
      #  delete.GP.SMED <- function(mod){}
      } else {
        message("Package not recognized Error # 1347344")
      }#;browser()
      if(length(X) != 0 & length(Z) != 0 & length(package) != 0) {
        init()
      }
    }, # end initialize
    init = function(X=NULL, Z=NULL, ...) {#browser()
      if (!is.null(X)) {X <<- X}
      if (!is.null(Z)) {Z <<- Z}
      if (length(.self$X) == 0 | length(.self$Z) == 0) {stop("X or Z not set")}
      n.at.last.update <<- nrow(.self$X)
      if (max(.self$Z) - min(.self$Z) < 1e-8) {warning("Z values are too close, adding noise"); .self$Z <<- .self$Z + rnorm(length(.self$Z), 0, 1e-6)}
      #mod <<- list(.init())
      .init(...=...)
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL, ...) {#browser()
      if (length(n.at.last.update) == 0) {
        init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
      } else {
        if (!is.null(Xall)) {X <<- Xall} else if (!is.null(Xnew)) {X <<- rbind(X, Xnew)}
        if (!is.null(Zall)) {Z <<- Zall} else if (!is.null(Znew)) {Z <<- c(Z, Znew)}
        .update(...=...)
      }
      n.at.last.update <<- nrow(X)
    }, # end update
    predict = function(XX, se.fit = FALSE, ...) {#browser()
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      .predict(XX, se.fit=se.fit, ...=...)
    },
    predict.se = function(XX, ...) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      .predict.se(XX, ...=...)
    },
    predict.var = function(XX, ...) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      .predict.var(XX, ...=...)
    },
    grad = function (XX) {#browser() # NUMERICAL GRAD IS OVER 10 TIMES SLOWER
      if (!is.matrix(XX)) {
        if (ncol(X) == 1) XX <- matrix(XX, ncol=1)
        else if (length(XX) == ncol(X)) XX <- matrix(XX, nrow=1)
        else stop('Predict input should be matrix')
      } else {
        if (ncol(XX) != ncol(X)) {stop("Wrong dimension input")}
      }
      grad.func <- function(xx) predict(xx)
      grad.apply.func <- function(xx) numDeriv::grad(grad.func, xx)
      grad1 <- apply(XX, 1, grad.apply.func)
      if (ncol(X) == 1) return(grad1)
      t(grad1)
    },
    grad_norm = function (XX) {#browser()
      grad1 <- grad(XX)
      if (!is.matrix(grad1)) return(abs(grad1))
      apply(grad1,1, function(xx) {sqrt(sum(xx^2))})
    },
    theta = function() {
      .theta()
    },
    nugget = function() {
      .nugget()
    },
    delete = function(...) {
      .delete(...=...)
    },
    finalize = function(...) {

    }
  )
)

