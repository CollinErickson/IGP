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
    mod = "list"
  ),
  methods = list(
    initialize = function(...) {#browser()
      callSuper(...)

      if (package == "GPfit") {
        .init <<- function() {GPfit::GP_fit(X, Z)}
        .update <<- function(){mod <<- list(GPfit::GP_fit(X, Z))}
        .predict <<- function(XX){GPfit::predict.GP(mod[[1]], XX)$Y_hat}
        .delete <<- function(){mod <<- list()}
      } else if (package=="laGP") {
        .init <<- function() {laGP::newGPsep(X=X, Z=Z, d=p, g=1e-8)}
        .update <<- function() {laGP::updateGPsep(gpsepi=mod, X=X, Z=Z);return(mod)}
        .predict <<- function(XX){laGP::predGPsep(mod, matrix(XX, 1, p))$mean}
        .delete <<- laGP::deleteGPsep
      } else if (package=="mlegp") {
        .init <<- function() {mlegp::mlegp(X=X, Z=Z, verbose=0)}
        .update <<- function() {mod <<- list(mlegp::mlegp(X=X, Z=Z, verbose=0))}
        .predict <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX)}
        .predict.se <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit}
        .predict.var <<- function(XX) {mlegp::predict.gp(object=mod[[1]], newData=XX, se.fit=T)$se.fit^2}
        .delete <<- function(){mod <<- list()}
      #} else if (GP.package=='exact') {
      #  predict.GP.SMED <- function(mod,xx) {f(xx)}
      #  init.GP.SMED <- function(X,Y) {}
      #  update.GP.SMED <- function(mod,X,Y) {}
      #  delete.GP.SMED <- function(mod){}
      } else {
        stop("No package specified Error # 579238572")
      }
      if(length(X) != 0 & length(Z) != 0) {
        init()
      }
    }, # end initialize
    init = function(X=NULL, Z=NULL) {#browser()
      if (!is.null(X)) {X <<- X}
      if (!is.null(Z)) {Z <<- Z}
      if (length(.self$X) == 0 | length(.self$Z) == 0) {stop("X or Z not set")}
      mod <<- list(.init())
      print('done')
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL) {
      if (!is.null(Xall)) {X <<- Xall}
      if (!is.null(Zall)) {Z <<- Zall}
      if (!is.null(Xnew)) {X <<- rbind(X, Xnew)}
      if (!is.null(Znew)) {X <<- c(Z, Znew)}
      .update()
    }, # end update
    predict = function(XX) {
      .predict(XX)
    },
    predict.se = function(XX) {
      .predict.se(XX)
    },
    predict.var = function(XX) {
      .predict.var(XX)
    },
    delete = function() {
      .delete()
    },
    finalize = function() {

    }
  )
  )

