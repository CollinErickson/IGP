UGP <- setRefClass("UGP",
  fields = list(
    X = "matrix",
    Z = "numeric",
    package = "character",
    init2 = "function",
    update2 = "function",
    predict2 = "function",
    predict.se2 = "function",
    predict.var2 = "function",
    delete2 = "function",
    mod = "list"
  ),
  methods = list(
    initialize = function(...) {#browser()
      callSuper(...)

      if (package == "GPfit") {
        init2 <<- GPfit::GP_fit
        update2 <<- function(){GPfit::GP_fit(X, Z)}
        predict2 <<- function(XX){GPfit::predict.GP(mod, matrix(XX, 1, p))$Y_hat}
        delete2 <<- function(){}
      } else if (package=="laGP") {
        init2 <<- function() {laGP::newGPsep(X=X, Z=Z, d=p, g=1e-8)}
        update2 <<- function() {laGP::updateGPsep(gpsepi=mod, X=X, Z=Z);return(mod)}
        predict2 <<- function(XX){laGP::predGPsep(mod, matrix(XX, 1, p))$mean}
        delete2 <<- laGP::deleteGPsep
      } else if (package=="mlegp") {
        init2 <<- function() {mlegp::mlegp(X=X, Z=Z, verbose=0)}
        update2 <<- function() {mod <<- list(mlegp::mlegp(X=X, Z=Z, verbose=0))}
        predict2 <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX)}
        predict.se2 <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit}
        predict.var2 <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit^2}
        delete2 <<- function(){rm(mod)}
      #} else if (GP.package=='exact') {
      #  predict.GP.SMED <- function(mod,xx) {f(xx)}
      #  init.GP.SMED <- function(X,Y) {}
      #  update.GP.SMED <- function(mod,X,Y) {}
      #  delete.GP.SMED <- function(mod){}
        if(length(X) != 0 & length(Z) != 0) {
          init()
        }
      } else {
        stop("No package specified Error # 579238572")
      }
    }, # end initialize
    init = function(X=NULL, Z=NULL) {#browser()
      if (!is.null(X)) {X <<- X}
      if (!is.null(Z)) {Z <<- Z}
      if (length(.self$X) == 0 | length(.self$Z) == 0) {stop("X or Z not set")}
      mod <<- list(init2())
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL) {
      if (!is.null(Xall)) {X <<- Xall}
      if (!is.null(Zall)) {Z <<- Zall}
      if (!is.null(Xnew)) {X <<- rbind(X, Xnew)}
      if (!is.null(Znew)) {X <<- c(Z, Znew)}
      update2()
    }, # end update
    predict = function(XX) {
      predict2(XX)
    },
    predict.se = function(XX) {
      predict.se2(XX)
    },
    predict.var = function(XX) {
      predict.var2(XX)
    },
    delete = function() {
      delete2()
    },
    finalize = function() {

    }
  )
  )

