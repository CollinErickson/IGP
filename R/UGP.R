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
    mod = "list",
    mod.extra = "list" # list to store additional data needed for model
  ),
  methods = list(
    initialize = function(...) {#browser()
      callSuper(...)

      if (length(package)==0) {
        #message("No package specified Error # 579238572")
      } else if (package == "GPfit") {
        .init <<- function() {GPfit::GP_fit(X, Z)}
        .update <<- function(){mod <<- list(GPfit::GP_fit(X, Z))}
        .predict <<- function(XX){GPfit::predict.GP(mod[[1]], XX)$Y_hat}
        .predict.se <<- function(XX) {sqrt(GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE)}
        .predict.var <<- function(XX) {GPfit::predict.GP(object=mod[[1]], xnew=XX, se.fit=T)$MSE}
        .delete <<- function(){mod <<- list()}
      } else if (package=="laGP") {
        .init <<- function() {
          da <- laGP::darg(list(mle=TRUE), X=X)
          ga <- laGP::garg(list(mle=TRUE), y=Z)
          mod.extra <<- list(da=da, ga=ga) # store extra data for update
          #laGP::newGPsep(X=X, Z=Z, d=p, g=1e-8)
          mod1 <- laGP::newGPsep(X=X, Z=Z, d=da$start, g=ga$start, dK = TRUE)
          laGP::jmleGPsep(gpsepi = mod1, drange=c(da$min, da$max),
                                 grange=c(ga$min, ga$max),
                                 dab=da$ab, gab=ga$ab, verb=1, maxit=1000)
          mod1
        }
        .update <<- function() {#browser()
          da <- mod.extra$da
          ga <- mod.extra$ga
          laGP::updateGPsep(gpsepi=mod[[1]], X=X, Z=Z)
          #mle <- laGP::jmleGPsep(gpsepi = mod[[1]], drange=c(da$min, da$max), grange=c(ga$min, ga$max), dab=da$ab, gab=ga$ab, verb=1)
          }
        .predict <<- function(XX){laGP::predGPsep(mod, XX, lite=TRUE)$mean}
        .delete <<- function() {laGP::deleteGPsep(mod[[1]]);mod <<- list()}
      } else if (package=="mlegp") {
        .init <<- function() {
          co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
          m
        }
        .update <<- function() {
          co <- capture.output(m <- mlegp::mlegp(X=X, Z=Z, verbose=0))
          mod <<- list(m)
        }
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
      mod <<- list(.init())
      print('done')
    }, # end init
    update = function(Xall=NULL, Zall=NULL, Xnew=NULL, Znew=NULL) {#browser()
      if (!is.null(Xall)) {X <<- Xall}
      if (!is.null(Zall)) {Z <<- Zall}
      if (!is.null(Xnew)) {X <<- rbind(X, Xnew)}
      if (!is.null(Znew)) {Z <<- c(Z, Znew)}
      .update()
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

