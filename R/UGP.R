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
    n.at.last.update = "numeric" # count how many in update, must be at end of X
  ),
  methods = list(
    initialize = function(...) {#browser()
      callSuper(...)

      if (length(package)==0) {
        #message("No package specified Error # 579238572")
      } else if (package == "GPfit") {
        .init <<- function() {mod <<- list(GPfit::GP_fit(X, Z))}
        .update <<- function(){mod <<- list(GPfit::GP_fit(X, Z))}
        .predict <<- function(XX){GPfit::predict.GP(mod[[1]], XX)$Y_hat}
        .predict.se <<- function(XX) {sqrt(GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE)}
        .predict.var <<- function(XX) {GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE}
        .delete <<- function(){mod <<- list()}
      } else if (package=="laGP") {
        .init <<- function() {
          da <- laGP::darg(list(mle=TRUE), X=X)
          ga <- laGP::garg(list(mle=TRUE), y=Z)
          mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=ga$start, dK = TRUE)
          #mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=1e-6, dK = TRUE)
          laGP::jmleGPsep(gpsepi = mod1, drange=c(da$min, da$max),
                                 grange=c(ga$min, ga$max),
                                 dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
          mod <<- list(mod1)
        }
        .update <<- function() {#browser()
          da <- laGP::darg(list(mle=TRUE), X=X)
          ga <- laGP::garg(list(mle=TRUE), y=Z)
          n.since.last.update = nrow(X) - n.at.last.update
          if (n.since.last.update < 1) {
            message("Can't update, no new X rows")
          } else {
            if (n.at.last.update < 10 || n.since.last.update > .25 * n.at.last.update) {
              # start over if too many
              .delete()
              .init()
            } else {
              laGP::updateGPsep(gpsepi=mod[[1]], X=X[-(1:n.at.last.update),], Z=Z[-(1:n.at.last.update)])
            }
          }
          laGP::jmleGPsep(gpsepi = mod[[1]], drange=c(da$min, da$max),
                          grange=c(ga$min, ga$max),
                          dab=da$ab, gab=ga$ab, verb=0, maxit=1000)
          }
        .predict <<- function(XX){laGP::predGPsep(mod, XX, lite=TRUE)$mean}
        .predict.se <<- function(XX) {sqrt(laGP::predGPsep(mod, XX, lite=TRUE)$s2)}
        .predict.var <<- function(XX) {laGP::predGPsep(mod, XX, lite=TRUE)$s2}
        .delete <<- function() {laGP::deleteGPsep(mod[[1]]);mod <<- list()}
      } else if (package=="mlegp") {
        .init <<- function() {
          co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
          mod <<- list(m)
        }
        .update <<- function() {
          co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
          mod <<- list(m)
        }
        .predict <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX)}
        .predict.se <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit}
        .predict.var <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit^2}
        .delete <<- function(){mod <<- list()}
      } else if (package == "GPy") {
        require("rPython")
        .init <<- function() {
          python.exec('import sys') # These first two lines need to go
          python.exec("sys.path.insert(0, '/Users/collin/anaconda/lib/python2.7/site-packages/')")
          python.exec('import numpy as np')
          python.exec('import GPy')

          python.assign("inputdim", ncol(X))
          python.assign("X1", (X))
          python.assign("y1", Z)
          python.exec('X =  np.matrix(X1)')
          python.exec('y = np.matrix(y1).reshape((-1,1))')
          python.exec("kernel = GPy.kern.RBF(input_dim=inputdim, variance=1., lengthscale=[1. for iii in range(inputdim)],ARD=True)")
          python.exec("gp = GPy.models.GPRegression(X,y,kernel,normalizer=True)")
          python.exec("gp.likelihood.variance = 1e-8")
          python.exec("gp.optimize(messages=False)")
          python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")

          mod <<- list("GPy model is in Python")
        }
        .update <<- function() {
          python.assign("X1", (X))
          python.assign("y1", Z)
          python.exec('X =  np.matrix(X1)')
          python.exec('y = np.matrix(y1).reshape((-1,1))')
          python.exec("gp.set_XY(X = X, Y = y)")
          python.exec("gp.optimize(messages=False)")
          python.exec("gp.optimize_restarts(num_restarts = 5,  verbose=False)")
        }
        .predict <<- function(XX) {
          python.assign("xp1", XX)
          python.exec("xp = np.asmatrix(xp1)")
          python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          unlist(python.get("y_pred.tolist()"))
        }
        .predict.se <<- function(XX) {
          python.assign("xp1", XX)
          python.exec("xp = np.asmatrix(xp1)")
          python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          unlist(python.get("np.sqrt(sigma2_pred).tolist()"))
        }
        .predict.var <<- function(XX) {
          python.assign("xp1", XX)
          python.exec("xp = np.asmatrix(xp1)")
          python.exec("y_pred, sigma2_pred = gp.predict(np.asarray(xp))")
          unlist(python.get("sigma2_pred.tolist()"))
        }
        .delete <<- function(){
          python.exec('X =  None')
          python.exec('y =  None')
          python.exec('xp =  None')
          python.exec('X1 =  None')
          python.exec('y1 =  None')
          python.exec('xp1 =  None')
          python.exec('y_pred =  None')
          python.exec('sigma2_pred =  None')
          python.exec('gp =  None')
          python.exec('kernel =  None')
          python.exec('inputdim =  None')
          mod <<- list()
          }
      #} else if (GP.package=='exact') {
      #  predict.GP.SMED <- function(mod,xx) {f(xx)}
      #  init.GP.SMED <- function(X,Y) {}
      #  update.GP.SMED <- function(mod,X,Y) {}
      #  delete.GP.SMED <- function(mod){}
      } else {
        message("Package not recognized Error # 1347344")
      }
      if(length(X) != 0 & length(Z) != 0 & length(package) != 0) {
        init()
      }
    }, # end initialize
    init = function(X=NULL, Z=NULL) {#browser()
      if (!is.null(X)) {X <<- X}
      if (!is.null(Z)) {Z <<- Z}
      if (length(.self$X) == 0 | length(.self$Z) == 0) {stop("X or Z not set")}
      n.at.last.update <<- nrow(.self$X)
      #mod <<- list(.init())
      .init()
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL) {#browser()
      if (length(n.at.last.update) == 0) {
        init(X = if(!is.null(Xall)) Xall else Xnew, Z = if (!is.null(Zall)) Zall else Znew)
      } else {
        if (!is.null(Xall)) {X <<- Xall} else if (!is.null(Xnew)) {X <<- rbind(X, Xnew)}
        if (!is.null(Zall)) {Z <<- Zall} else if (!is.null(Znew)) {Z <<- c(Z, Znew)}
        .update()
      }
      n.at.last.update <<- nrow(X)
    }, # end update
    predict = function(XX) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      .predict(XX)
    },
    predict.se = function(XX) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      .predict.se(XX)
    },
    predict.var = function(XX) {
      if(!is.matrix(XX)) XX <- matrix(XX,nrow=1)
      .predict.var(XX)
    },
    delete = function() {
      .delete()
    },
    finalize = function() {

    }
  )
  )

